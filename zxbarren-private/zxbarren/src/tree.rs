use num::Rational;
use quizx::circuit::*;

fn add_ry(c: &mut Circuit, i: usize, param: &String) {
    c.add_gate_with_phase("rz", vec![i], Rational::new(-1, 2));
    c.add_gate_with_param("rx", vec![i], param);
    c.add_gate_with_phase("rz", vec![i], Rational::new(1, 2));
}

pub fn tree(n: usize) -> (Circuit, Vec<String>) {
    let mut c = Circuit::new(n);
    let mut params: Vec<String> = Vec::new();

    let mut qubits: Vec<usize> = (0..n).collect();
    let mut layer = 0;
    while qubits.len() > 0 {
        for &i in &qubits {
            let p = format!("alpha_{layer}_{i}");
            add_ry(&mut c, i, &p);
            params.push(p);
        }

        let mut new_qubits = Vec::new();
        let mut i = 0;
        while i < qubits.len() -1 {
            c.add_gate("cx", vec![qubits[i], qubits[i+1]]);
            new_qubits.push(qubits[i+1]);
            i += 2;
        }
        qubits = new_qubits;
        layer += 1;
    }

    (c, params)
}
