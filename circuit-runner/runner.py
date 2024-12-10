# Benchmark runtime: 23313.70 seconds = 6.48 hours

import my_pyzx as zx
from my_pyzx.basicrules import fuse
from my_pyzx.graph.base import BaseGraph
from fractions import Fraction
import numpy as np
import sympy as sp
from collections import defaultdict
import matplotlib.pyplot as plt
import random
import copy
import sys
import subprocess
import multiprocessing
import threading
import time

np.set_printoptions(threshold=np.inf) # print the whole matrix

B = zx.VertexType.BOUNDARY
Z = zx.VertexType.Z
X = zx.VertexType.X
H_BOX = zx.VertexType.H_BOX
W_INP = zx.VertexType.W_INPUT
W_OUT = zx.VertexType.W_OUTPUT
SE = zx.EdgeType.SIMPLE
HE = zx.EdgeType.HADAMARD

draw_scale_default = 30

def nice_print(str, conditional=True):
    if conditional:
        print(str)

# like zx.draw but with my favorite settings already set
def nice_draw(g, conditional=True, labels=True, scale=draw_scale_default):
    if not conditional:
        return
    
    exact_string = str(g.scalar).split(" = ")[1]
    numeric_string = str(g.scalar.to_number())
    print("\t" + exact_string + " ≈ " + numeric_string)

    zx.draw(g, labels=labels, scale=scale)

def to_pure_matrix(matrix):
    # get rid of floating point error
    # change -0.0 to 0.0
    return np.around(matrix, 10) + 0.

class RandQCirc:
    def __init__(self, qbit_cnt, not_cnt, allow_not_in_targ_range, cnot_cnt, allow_cnot_in_targ_range, tof_cnt, tof_targ_qbit_range, allow_ctrl_in_targ_range, seed=False):        
        self.qbit_cnt = qbit_cnt
        self.not_cnt = not_cnt
        self.cnot_cnt = cnot_cnt
        self.tof_cnt = tof_cnt
        self.allow_not_in_targ_range = allow_not_in_targ_range
        self.allow_cnot_in_targ_range = allow_cnot_in_targ_range
        self.tof_targ_qbit_range = tof_targ_qbit_range
        self.allow_ctrl_in_targ_range = allow_ctrl_in_targ_range

        self.string_circ_repr = [] # contains lists which represents column by column of the generated quantum circuit
        self.g = zx.Graph()

        for i in range(self.qbit_cnt):
            b1 = self.g.add_vertex(B, qubit=i, row=0)
            b2 = self.g.add_vertex(B, qubit=i, row=1)
            self.g.add_edge((b1, b2))
        self.g.auto_detect_io()

        if seed:
            random.seed(seed)

    # this creates the self.string_circ_repr which later will be used to convert it into a diagram
    def generate_string_circuit(self):
        # returns a list with self.qbit_cnt length, where the region enclosed formed by indices from self.tof_targ_qbit_range holds
        # the string "targ" in a random place, the other places are "". in the other region, there is a random amount of
        # strings "ctrl", however at least two and the other places are also ""
        def rand_tof_qbit_order() -> list:
            tof_qbit_order = [""]*self.qbit_cnt
            tof_qbit = random.choice(self.tof_targ_qbit_range)
            tof_qbit_order[tof_qbit] = "targ"

            if self.allow_ctrl_in_targ_range:
                pot_ctrl_indices = [i for i in range(self.qbit_cnt) if i != tof_qbit]
            else:
                pot_ctrl_indices = [i for i in range(self.qbit_cnt) if not i in self.tof_targ_qbit_range]
            
            num_of_ctrls = random.randint(2, len(pot_ctrl_indices))
            
            actual_ctrl_indices = random.sample(pot_ctrl_indices, num_of_ctrls)
            for i in actual_ctrl_indices:
                tof_qbit_order[i] = "ctrl"
            
            return tof_qbit_order
        
        # returns a list with self.qbit_cnt length, where one of the elements is "Z" and one of them "X", the rest ""
        def rand_cnot_qbit_order() -> list:
            cnot_qbit_order = [""]*self.qbit_cnt

            if self.allow_cnot_in_targ_range:
                cnot_qbit_order[random.choice(range(self.qbit_cnt))] = "Z"
                cnot_qbit_order[random.choice([i for i in range(self.qbit_cnt) if cnot_qbit_order[i] != "Z"])] = "X"
            else:
                cnot_qbit_order[random.choice([i for i in range(self.qbit_cnt) if i not in self.tof_targ_qbit_range])] = "Z"
                cnot_qbit_order[random.choice([i for i in range(self.qbit_cnt) if cnot_qbit_order[i] != "Z" and i not in self.tof_targ_qbit_range])] = "X"

            return cnot_qbit_order
        
        # returns a list with self.qbit_cnt length, where one of the elements is "Z" and one of them "X", the rest ""
        def rand_not_qbit_order() -> list:
            not_qbit_order = [""]*self.qbit_cnt

            if self.allow_not_in_targ_range:
                not_qbit_order[random.choice(range(self.qbit_cnt))] = "NOT"
            else:
                not_qbit_order[random.choice([i for i in range(self.qbit_cnt) if i not in self.tof_targ_qbit_range])] = "NOT"

            return not_qbit_order

        # creates a list of columns (column represented by a string like ["", "ctr", "ctrl", "", "targ"])
        gate_order = ["not"] * self.not_cnt + ["cnot"] * self.cnot_cnt + ["tof"] * self.tof_cnt
        random.shuffle(gate_order)

        for gate in gate_order:
            if gate == "not":
                self.string_circ_repr.append(rand_not_qbit_order())     
            elif gate == "cnot":
                self.string_circ_repr.append(rand_cnot_qbit_order())          
            elif gate == "tof":
                self.string_circ_repr.append(rand_tof_qbit_order())
    
    # prints a matrix format of self.string_circ_repr
    def print_string_circ_repr(self):
        row_list = list(zip(*self.string_circ_repr))

        for row in row_list:
            for entry in row:
                entry = entry if entry != "" else "#"
                print(entry, end="\t")
            print()

    # adds a NOT to self.g
    def add_not(self, qbit):
        cur = zx.Graph()

        for i in range(self.qbit_cnt):
            b1 = cur.add_vertex(B, qubit=i, row=0)

            if i == qbit:
                x = cur.add_vertex(X, qubit=i, row=1, phase=1)
                b2 = cur.add_vertex(B, qubit=i, row=2)
                cur.add_edge((b1, x))
                cur.add_edge((x, b2))
            else:
                b2 = cur.add_vertex(B, qubit=i, row=2)
                cur.add_edge((b1, b2))

        cur.auto_detect_io()

        self.g = self.g + cur

    # adds a CNOT to self.g
    def add_cnot(self, z_qbit, x_qbit):
        cur = zx.Graph()

        for i in range(self.qbit_cnt):
            b1 = cur.add_vertex(B, qubit=i, row=0)

            if i == z_qbit:
                z = cur.add_vertex(Z, qubit=i, row=1)
                b2 = cur.add_vertex(B, qubit=i, row=2)
                cur.add_edge((b1, z))
                cur.add_edge((z, b2))
            elif i == x_qbit:
                x = cur.add_vertex(X, qubit=i, row=1)
                b2 = cur.add_vertex(B, qubit=i, row=2)
                cur.add_edge((b1, x))
                cur.add_edge((x, b2))
            else:
                b2 = cur.add_vertex(B, qubit=i, row=2)
                cur.add_edge((b1, b2))

        cur.add_edge((z, x))
        cur.auto_detect_io()

        self.g = self.g + cur

    # adds a multi-controlled Toffoli gate to self.g
    def add_toffoli(self, ctrl_qbits, targ_qbit):
        cur = zx.Graph()

        cur = zx.Graph()

        z_verts = [] # notice this will be sorted from smallest to biggest
        x_vert = None

        for i in range(self.qbit_cnt):
            if i in ctrl_qbits:
                start = cur.add_vertex(B, qubit=i, row=0)
                z = cur.add_vertex(Z, qubit=i, row=1)
                cur.set_vdata(z, "label", "pot_master_node")
                z_verts.append(z)
                end = cur.add_vertex(B, qubit=i, row=3)
                cur.add_edge((start, z))
                cur.add_edge((z, end))
            elif i == targ_qbit:
                start = cur.add_vertex(B, qubit=i, row=0)
                x = cur.add_vertex(X, qubit=i, row=1)
                x_vert = x
                end = cur.add_vertex(B, qubit=i, row=3)
                cur.add_edge((start, x))
                cur.add_edge((x, end))
            else:
                start = cur.add_vertex(B, qubit=i, row=0)
                end = cur.add_vertex(B, qubit=i, row=3)
                cur.add_edge((start, end))
        
        targ_shift = cur.add_vertex(Z, qubit=-0.5, row=3, phase=Fraction(1,1))
        z_hbox = cur.add_vertex(Z, qubit=-1-0.2, row=3)
        cur.set_vdata(z_hbox, "label", "star_node_z")
        hbox = cur.add_vertex(H_BOX, qubit=-1, row=3)
        cur.scalar.add_power(-2)
        cur.add_edge((hbox, z_hbox))
        b_shift = cur.add_vertex(B, qubit=-1, row=3.5)
        cur.add_edge((x_vert, targ_shift))
        cur.add_edge((targ_shift, hbox))
        cur.add_edge((hbox, b_shift))

        for i, z_vert in enumerate(z_verts):
            x_pi = cur.add_vertex(X, qubit=-len(z_verts)+i-1, row=3, phase=1)
            z_hbox = cur.add_vertex(Z, qubit=-len(z_verts)+i-1-0.2, row=3.5)
            cur.set_vdata(z_hbox, "label", "star_node_z")
            hbox = cur.add_vertex(H_BOX, qubit=-len(z_verts)+i-1, row=3.5, phase=Fraction(1,1))
            cur.scalar.add_power(-2)
            cur.add_edge((hbox, z_hbox))
            b_shift = cur.add_vertex(B, qubit=-len(z_verts)+i-1, row=4)
            cur.add_edge((z_vert, x_pi))
            cur.add_edge((x_pi, hbox))
            cur.add_edge((hbox, b_shift))

        end_part = zx.Graph()
        mid_idx = len(ctrl_qbits + [targ_qbit])/2
        mid_idx = round(mid_idx * 2) / 2 # round to nearest 0.5
        if int(mid_idx) == mid_idx: # add 0.5 if ending in .0
            mid_idx += 0.5
        master = end_part.add_vertex(Z, qubit=-mid_idx, row=1, phase=Fraction(1,1))
        end_part.set_vdata(master, "label", "pot_master_node")
        final_not = end_part.add_vertex(X, qubit=-1, row=1, phase=Fraction(1,1))

        for i in range(len(z_verts)+1):
            bound = end_part.add_vertex(B, qubit=-(len(z_verts)+1)+i, row=0)
            if i == len(z_verts):
                end_part.add_edge((bound, final_not))
                end_part.add_edge((final_not, master))
            else:
                end_part.add_edge((bound, master))
        
        for i in range(self.qbit_cnt):
            end_part = end_part @ zx.generate.identity(1)
        
        cur.auto_detect_io()
        end_part.auto_detect_io()
        cur = cur + end_part
        cur.auto_detect_io()

        self.g = self.g + cur

    # generates self.g from the already generated self.string_circ_repr
    def generate_diagram(self) -> BaseGraph:
        for column in self.string_circ_repr:
            if "NOT" in column:
                self.add_not(column.index("NOT"))
            elif "Z" in column:
                self.add_cnot(column.index("Z"), column.index("X"))
            elif "ctrl" in column:
                self.add_toffoli([i for i in range(len(column)) if column[i] == "ctrl"], column.index("targ"))

def convert_had_to_pi_ov_3_placeholder(g):
    for v in g.clone().vertices():
        if "label" in g.vdata_keys(v):
            if "star_node_z" == g.vdata(v, "label"):
                n = list(g.neighbors(v))[0]
                g.remove_vertex(v)
                
                n_qubit = g.qubit(n)
                n_row = g.row(n)

                n_n = list(g.neighbors(n))
                g.remove_vertex(n)

                placeholder = g.add_vertex(Z, qubit=n_qubit, row=n_row, phase=Fraction(1,3))

                g.add_edge((n_n[0], placeholder))
                g.add_edge((n_n[1], placeholder))

def to_dot(g) -> str:
    dot = "graph {\n"
    
    for v in g.vertices():
        t = g.type(v)
        p = g.phase(v)
        
        color = {
            B: "black",
            Z: "green",
            X: "red"
        }.get(t, "black")
        
        label = ""
        if v in g.inputs():
            label = f"{v}:i"
        elif v in g.outputs():
            label = f"{v}:o"
        elif p != 0:
            label = f"{v}:{p}"
        else:
            label = f"{v}"
        
        dot += f"  {v} [color={color}, label=\"{label}\""
        
        q = g.qubit(v)
        r = g.row(v)
        if q != 0 or r != 0:
            dot += f", pos=\"{q},{r}!\""
        
        dot += "]\n"
    
    dot += "\n"
    
    for src, trg in g.edges():
        e = g.edge(src, trg)
        et = g.edge_type(e)

        dot += f"  {src} -- {trg}"
        if et == HE:
            dot += " [color=blue]"
        # if ty == EType.T:
        #     dot += " [color=orange]"
        dot += "\n"
    
    dot += "}\n"
    
    return dot

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

timeout = 180

def thread1_quizx_with(qbit, not_cnt, cnot_cnt, tof_cnt, my_seed, result_dict):
    command = ["./target/release/main", "WITH_SD", "benchmark_dot_repr/{}_{}_{}_{}_{}.txt".format(qbit, not_cnt, cnot_cnt, tof_cnt, my_seed)]
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
        result_dict["result_quizx_with"] = result.stdout
    except subprocess.TimeoutExpired:
        result_dict["result_quizx_with"] = "QUIZX_WITH NUM OF TERMS: NONE\nQUIZX_WITH RUNTIME: STOPPED\n\n" # oopsie, bad newlines

def thread2_quizx_without(qbit, not_cnt, cnot_cnt, tof_cnt, my_seed, result_dict):
    command = ["./target/release/main", "WITHOUT_SD","benchmark_dot_repr/{}_{}_{}_{}_{}.txt".format(qbit, not_cnt, cnot_cnt, tof_cnt, my_seed)]
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
        result_dict["result_quizx_without"] = result.stdout
    except subprocess.TimeoutExpired:
        result_dict["result_quizx_without"] = "QUIZX_WITHOUT NUM OF TERMS: NONE\nQUIZX_WITHOUT RUNTIME: STOPPED\n\n" # oopsie, bad newlines


def program(args):
    qbit, not_cnt, cnot_cnt, tof_cnt, my_seed = args

    randqcirc = RandQCirc(qbit_cnt=qbit, not_cnt=not_cnt, allow_not_in_targ_range=False, cnot_cnt=cnot_cnt, tof_cnt=tof_cnt, allow_cnot_in_targ_range=False, tof_targ_qbit_range=range(qbit-int(0.15*qbit),qbit), allow_ctrl_in_targ_range=False, seed=my_seed)
    randqcirc.generate_string_circuit()
    randqcirc.generate_diagram()
    randqcirc.g.apply_state("+"*(qbit-int(0.15*qbit)) + "0"*int(0.15*qbit))
    randqcirc.g.apply_effect("+"*(qbit-int(0.15*qbit)) + "0"*int(0.15*qbit))
    g = randqcirc.g

    g_quizx = g.clone()
    convert_had_to_pi_ov_3_placeholder(g_quizx)
    dot_string = to_dot(g_quizx).rstrip()
    with open("benchmark_dot_repr/{}_{}_{}_{}_{}.txt".format(qbit, not_cnt, cnot_cnt, tof_cnt, my_seed), "w") as file:
        file.write(dot_string)

    result_dict = {}

    thread1 = threading.Thread(target=thread1_quizx_with, args=(qbit, not_cnt, cnot_cnt, tof_cnt, my_seed, result_dict))
    thread2 = threading.Thread(target=thread2_quizx_without, args=(qbit, not_cnt, cnot_cnt, tof_cnt, my_seed, result_dict))
    
    thread1.start()
    thread2.start()
    
    # wait for threads to finish
    thread1.join()
    thread2.join(timeout=timeout)

    with open("benchmark_output/{}_{}_{}_{}_{}.txt".format(qbit, not_cnt, cnot_cnt, tof_cnt, my_seed), "w") as file:
        file.write("{}_{}_{}_{}_{}\n".format(qbit, not_cnt, cnot_cnt, tof_cnt, my_seed) + result_dict["result_quizx_with"] + result_dict["result_quizx_without"])

def main():
    input("Press Enter to start...")
    
    t0 = time.time()

    num_processes = 12  # number of processes to run concurrently = number of cores

    tasks = []
    for qbit in [20, 50, 80]:
        for not_cnt in [0, 20, 80, 240]:
            for cnot_cnt in [0, 20, 80, 240]:
                for tof_cnt in [20, 40, 60, 80]:
                    for sample_i in range(50):
                        my_seed = random.randint(1, 100000000)
                        tasks.append([qbit, not_cnt, cnot_cnt, tof_cnt, my_seed])

    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(program, tasks)

    t1 = time.time()

    print("Benchmark runtime:", t1-t0)

if __name__ == "__main__":
    main()