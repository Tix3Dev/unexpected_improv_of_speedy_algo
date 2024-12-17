use quizx::circuit::*;
use num::rational::Rational;


pub fn circ1(n: usize, l: usize) -> (Circuit, Vec<String>) {
    let mut c = Circuit::new(n);
    let mut params: Vec<String> = Vec::new();

    for k in 0..l {
        for i in 0..n {
            let p1 = format!("alpha_{k}_{i}");
            let p2 = format!("beta_{k}_{i}");
            c.add_gate_with_param("rx", vec![i], &p1);
            c.add_gate_with_param("rz", vec![i], &p2);
            params.push(p1);
            params.push(p2);
        }
    }

    (c, params)
}


pub fn circ2(n: usize, l: usize) -> (Circuit, Vec<String>) {
    let mut c = Circuit::new(n);
    let mut params: Vec<String> = Vec::new();

    for k in 0..l {
        for i in 0..n {
            let p1 = format!("alpha_{k}_{i}");
            let p2 = format!("beta_{k}_{i}");
            c.add_gate_with_param("rx", vec![i], &p1);
            c.add_gate_with_param("rz", vec![i], &p2);
            params.push(p1);
            params.push(p2);
        }
        let mut i = n-1;
        while i > 0 {
            c.add_gate("cx", vec![i, i-1]);
            i -= 1
        }
    }

    (c, params)
}


pub fn circ9(n: usize, l: usize) -> (Circuit, Vec<String>) {
    let mut c = Circuit::new(n);
    let mut params: Vec<String> = Vec::new();

    for k in 0..l {
        for i in 0..n {
            c.add_gate("h", vec![i]);
        }
        for i in 0..n-1 {
            c.add_gate("cz", vec![i, i+1]);
        }
        for i in 0..n {
            let p = format!{"alpha_{k}_{i}"};
            c.add_gate_with_param("rx", vec![i], &p);
            params.push(p);
        }
    }
    
    (c, params)
}

fn add_ry(c: &mut Circuit, i: usize, param: &String) {
    c.add_gate_with_phase("rz", vec![i], Rational::new(-1, 2));
    c.add_gate_with_param("rx", vec![i], param);
    c.add_gate_with_phase("rz", vec![i], Rational::new(1, 2));
}


pub fn circ10(n: usize, l: usize) -> (Circuit, Vec<String>) {
    let mut c = Circuit::new(n);
    let mut params: Vec<String> = Vec::new();

    for i in 0..n {
        let p = format!{"alpha_{i}"};
        add_ry(&mut c, i, &p);
        params.push(p);
    }

    for k in 0..l {
        for i in 0..n-1 {
            c.add_gate("cz", vec![i, i+1]);
        }
        c.add_gate("cz", vec![0, n-1]);
        for i in 0..n {
            let p = format!{"beta_{k}_{i}"};
            add_ry(&mut c, i, &p);
            params.push(p);
        }
    }
    
    (c, params)
}


pub fn circ11(n: usize, l: usize) -> (Circuit, Vec<String>) {
    let mut c = Circuit::new(n);
    let mut params: Vec<String> = Vec::new();

    for k in 0..l {
        for i in 0..n {
            let p1 = format!{"alpha_{k}_{i}"};
            let p2 = format!{"beta_{k}_{i}"};
            add_ry(&mut c, i, &p1);
            c.add_gate_with_param("rz", vec![i], &p2);
            params.push(p1);
            params.push(p2);
        }

        let mut i = 0;
        while i < n-1 {
            c.add_gate("cx", vec![i, i+1]);
            i += 2;
        }

        let mut i = 1;
        while i < n {
            let p1 = format!{"gamma_{k}_{i}"};
            let p2 = format!{"delta_{k}_{i}"};
            c.add_gate_with_param("ry", vec![i], &p1);
            c.add_gate_with_param("rz", vec![i], &p2);
            params.push(p1);
            params.push(p2);
            if i+1 < n {
                let p3 = format!{"epsilon_{k}_{i}"};
                let p4 = format!{"phi_{k}_{i}"};
                add_ry(&mut c, i+1, &p3);
                c.add_gate_with_param("rz", vec![i+1], &p4);
                c.add_gate("cx", vec![i, i+1]);
                params.push(p3);
                params.push(p4);
            }
            i += 3;
        }
    }
    
    (c, params)
}


pub fn circ12(n: usize, l: usize) -> (Circuit, Vec<String>) {
    let mut c = Circuit::new(n);
    let mut params: Vec<String> = Vec::new();

    for k in 0..l {
        for i in 0..n {
            let p1 = format!{"alpha_{k}_{i}"};
            let p2 = format!{"beta_{k}_{i}"};
            add_ry(&mut c, i, &p1);
            c.add_gate_with_param("rz", vec![i], &p2);
            params.push(p1);
            params.push(p2);
        }

        let mut i = 0;
        while i < n-1 {
            c.add_gate("cz", vec![i, i+1]);
            i += 2;
        }

        let mut i = 1;
        while i < n {
            let p1 = format!{"gamma_{k}_{i}"};
            let p2 = format!{"delta_{k}_{i}"};
            add_ry(&mut c, i, &p1);
            c.add_gate_with_param("rz", vec![i], &p2);
            params.push(p1);
            params.push(p2);
            if i+1 < n {
                let p3 = format!{"epsilon_{k}_{i}"};
                let p4 = format!{"phi_{k}_{i}"};
                add_ry(&mut c, i+1, &p3);
                c.add_gate_with_param("rz", vec![i+1], &p4);
                c.add_gate("cx", vec![i, i+1]);
                params.push(p3);
                params.push(p4);
            }
            i += 3;
        }
    }
    
    (c, params)
}


pub fn circ15(n: usize, l: usize) -> (Circuit, Vec<String>) {
    let mut c = Circuit::new(n);
    let mut params: Vec<String> = Vec::new();

    for k in 0..l {
        for i in 0..n {
            let p = format!{"alpha_{k}_{i}"};
            add_ry(&mut c, i, &p);
            params.push(p);
        }

        c.add_gate("cx", vec![0, n-1]);

        let mut i = n-1;
        while i > 0 {
            c.add_gate("cx", vec![i, i-1]);
            i -= 1;
        }

        for i in 0..n {
            let p = format!{"beta_{k}_{i}"};
            add_ry(&mut c, i, &p);
            params.push(p);
        }

        c.add_gate("cx", vec![n-2, n-1]);
        c.add_gate("cx", vec![n-1, 0]);
        for i in 1..n-1 {
            c.add_gate("cx", vec![i-1, i]);
        }
    }
    
    (c, params)
}

