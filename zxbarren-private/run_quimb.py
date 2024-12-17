import sys
from time import time

from pyzx import Circuit
from pyzx.circuit.gates import ZPhase, Z, NOT, CNOT, HAD, CZ
from pyzx.simplify import to_gh
from pyzx.quimb import to_quimb_tensor
from sympy import Symbol, Integer
from sympy import default_sort_key

def sim1(n, l):
    circ = Circuit(qubit_amount=n)

    for k in range(l):
        for i in range(n):
            p1 = Symbol(f'alpha_{k}_{i}')
            p2 = Symbol(f'beta_{k}_{i}')
            circ.add_gate(HAD(i))
            circ.add_gate(ZPhase(i, p1))
            circ.add_gate(HAD(i))
            circ.add_gate(ZPhase(i, p2))
    return circ

def sim2(n, l):
    circ = Circuit(qubit_amount=n)

    for k in range(l):
        for i in range(n):
            p1 = Symbol(f'alpha_{k}_{i}')
            p2 = Symbol(f'beta_{k}_{i}')
            circ.add_gate(HAD(i))
            circ.add_gate(ZPhase(i, p1))
            circ.add_gate(HAD(i))
            circ.add_gate(ZPhase(i, p2))
        for i in range(n - 1):
            circ.add_gate(CNOT(i, i + 1))
    return circ

def sim9(n, l):
    circ = Circuit(qubit_amount=n)

    for k in range(l):
        for i in range(n):
            circ.add_gate(HAD(i))
        for i in range(n - 1):
            circ.add_gate(CZ(i, i + 1))
        for i in range(n):
            circ.add_gate(HAD(i))
            circ.add_gate(ZPhase(i, Symbol(f'alpha_{k}_{i}')))
            circ.add_gate(HAD(i))
    return circ

def add_ry(circ, i, param):
    circ.add_gate(ZPhase(i, -1/Integer(2)))
    circ.add_gate(ZPhase(i, param))
    circ.add_gate(ZPhase(i, 1/Integer(2)))

def sim10(n, l):
    circ = Circuit(qubit_amount=n)
    
    for i in range(n):
        add_ry(circ, i, Symbol(f'alpha_{i}'))
    for k in range(l):
        for i in range(n-1):
            circ.add_gate(CZ(i, i + 1))
        circ.add_gate(CZ(0, n - 1))
        for i in range(n - 1):
            circ.add_gate(CZ(i, i + 1))
        for i in range(n):
            add_ry(circ, i, Symbol(f'beta_{k}_{i}'))
    return circ

def sim11(n, l):
    circ = Circuit(qubit_amount=n)

    for k in range(l):
        for i in range(n):
            add_ry(circ, i, Symbol(f'alpha_{k}_{i}'))
            circ.add_gate(ZPhase(i, Symbol(f'beta_{k}_{i}')))
            
        i = 0
        while i < n-1:
            circ.add_gate(CNOT(i, i+1))
            i += 2
        
        i = 1
        while i < n:
            add_ry(circ, i, Symbol(f'gamma_{k}_{i}'))
            circ.add_gate(ZPhase(i, Symbol(f'delta_{k}_{i}')))
            if i+1 < n:
                add_ry(circ, i+1, Symbol(f'epsilon_{k}_{i}'))
                circ.add_gate(ZPhase(i+1, Symbol(f'phi_{k}_{i}')))
                circ.add_gate(CNOT(i, i+1))
            i += 3
    return circ

def sim12(n, l):
    circ = Circuit(qubit_amount=n)

    for k in range(l):
        for i in range(n):
            add_ry(circ, i, Symbol(f'alpha_{k}_{i}'))
            circ.add_gate(ZPhase(i, Symbol(f'beta_{k}_{i}')))
            
        i = 0
        while i < n-1:
            circ.add_gate(CZ(i, i+1))
            i += 2
        
        i = 1
        while i < n:
            add_ry(circ, i, Symbol(f'gamma_{k}_{i}'))
            circ.add_gate(ZPhase(i, Symbol(f'delta_{k}_{i}')))
            if i+1 < n:
                add_ry(circ, i+1, Symbol(f'epsilon_{k}_{i}'))
                circ.add_gate(ZPhase(i+1, Symbol(f'phi_{k}_{i}')))
                circ.add_gate(CZ(i, i+1))
            i += 3
    return circ

def sim15(n, l):
    circ = Circuit(qubit_amount=n)

    for k in range(l):
        for i in range(n):
            add_ry(circ, i, Symbol(f'alpha_{k}_{i}'))
        circ.add_gate(CNOT(0, n-1))
            
        i = n-1
        while i > 0:
            circ.add_gate(CNOT(i, i-1))
            i -= 2
            
        for i in range(n):
            add_ry(circ, i, Symbol(f'beta_{k}_{i}'))
        circ.add_gate(CNOT(n-2, n-1))
        circ.add_gate(CNOT(n-1, 0))
        
        for i in range(1,n-1):
            circ.add_gate(CNOT(i-1, i))
    return circ

def hamiltonian(paulis):
    circ = Circuit(qubit_amount=len(paulis))
    for i, p in enumerate(paulis):
        if p == 'Z':
            circ.add_gate(Z(i))
        elif p == 'X':
            circ.add_gate(NOT(i))
        elif p == 'Y':
            circ.add_gate(Z(i))
            circ.add_gate(NOT(i))
        elif p != 'I':
            raise ValueError(f'Invalid Pauli {p}')
    return circ

def get_params(g):
    exprs = [e for e in g.phases().values() if hasattr(e, 'free_symbols')]
    return sorted(
        set.union(set(), *[e.free_symbols for e in exprs]),
        key=default_sort_key)

def variance(ansatz, ham, param_index):
    exp = (ansatz + ham + ansatz.adjoint()).to_graph()
    Z_TYPE = 1
    X_TYPE = 2
    for i in exp.inputs():
        exp.set_type(i, X_TYPE)
    for i in exp.outputs():
        exp.set_type(i, X_TYPE)
    var = exp.copy()
    var.merge(exp)
    # to_gh(var)
    phase_vs = sorted([
        (k, v)
        for k, v in var.phases().items()
        if hasattr(v, 'free_symbols')
    ], key=lambda x: x[0])
    syms = get_params(var)
    for i, sym in enumerate(syms):
        v0, v1, v2, v3 = [k for k, v in phase_vs if sym in v.free_symbols]
        if i == param_index:
            for v in [v0, v1, v2, v3]:
                var.set_phase(v, 0)
                var.set_type(v, Z_TYPE)
            for x, y in [(v0, v1), (v0, v2), (v1, v3)]:
                p = var.add_vertex(X_TYPE, phase=1)
                var.add_edge((p, x))
                var.add_edge((p, y))
        else:
            for v in [v0, v1, v2, v3]:
                var.set_phase(v, 0)
            nv2 = var.add_vertex(X_TYPE, phase=0)
            nv3 = var.add_vertex(X_TYPE, phase=0)
            var.add_edge((v2, nv2))
            var.add_edge((v3, nv3))

            var.add_edge((v0, nv2))
            var.add_edge((v1, nv3))
            v4 = var.add_vertex(X_TYPE, phase=0)
            var.add_edge((v0, v4))
            var.add_edge((v1, v4))
            v5 = var.add_vertex(Z_TYPE, phase=0)
            var.add_edge((nv2, v5))
            var.add_edge((nv3, v5))
            t1 = var.add_vertex(X_TYPE, phase=0)
            t2 = var.add_vertex(Z_TYPE, phase=-1/4)
            t3 = var.add_vertex(Z_TYPE, phase=1/4)
            t4 = var.add_vertex(X_TYPE, phase=0)
            t5 = var.add_vertex(Z_TYPE, phase=-1/4)
            t6 = var.add_vertex(X_TYPE, phase=0)
            t7 = var.add_vertex(Z_TYPE, phase=1/4)
            t8 = var.add_vertex(Z_TYPE, phase=0)
            var.add_edge((v4, t1))
            var.add_edge((t1, t2))
            var.add_edge((t1, t3))
            var.add_edge((t2, t4))
            var.add_edge((t4, t5))
            var.add_edge((t3, t6))
            var.add_edge((t6, t7))
            var.add_edge((t4, t8))
            var.add_edge((t6, t8))
            var.add_edge((v5, t8))
            var.scalar.power2 -= 2  # divide by 2
    # normalise all of the X spiders by 2^((m + n - 2) / 2)
    for v in var.vertices():
        if var.type(v) == X_TYPE:
            n_neighbours = len(list(var.neighbors(v)))
            var.scalar.power2 += n_neighbours - 2
    return var



def run(name, n, l, backend="numpy", parallel=False):
    if name == 'sim2':
        ansatz = sim2(n, l)
    if name == 'sim9':
        ansatz = sim9(n, l)
    if name == 'sim10':
        ansatz = sim10(n, l)
    if name == 'sim11':
        ansatz = sim11(n, l)
    if name == 'sim12':
        ansatz = sim12(n, l)
    if name == 'sim15':
        ansatz = sim15(n, l)

    g = variance(ansatz, hamiltonian('Z' * n), 0)
    t = to_quimb_tensor(g)

    ts = time()

    t.contract(output_inds=[], backend=backend)

    contract_time = time() - ts
    print(contract_time)


if __name__ == "__main__":
    args = sys.argv[1:]
    run(args[0], int(args[1]), int(args[2]), 'jax' if 'jax' in args[3] else 'numpy', bool(args[4]))

