use zxbarren::variance::*;
use zxbarren::sim;


fn run_test(circ: &str, n: usize, l: usize, param: usize, hamiltonian: &str, expected: &str) {
    let (c, params) = match circ {
        "sim1" => sim::circ1(n, l),
        "sim2" => sim::circ2(n, l),
        "sim9" => sim::circ9(n, l),
        "sim10" => sim::circ10(n, l),
        "sim11" => sim::circ11(n, l),
        "sim12" => sim::circ12(n, l),
        "sim15" => sim::circ15(n, l),
        _ => panic!("Circuit not found: {circ}")
    };
    
    let var = compute_variance(&c, hamiltonian, &params[param]);
    assert_eq!(format!("{var}"), *expected);
}



#[test]
fn sim1_1() {
    run_test("sim1", 4, 1, 0, "ZZZZ", "2^-4 * (1)");
}

#[test]
fn sim1_2() {
    run_test("sim1", 4, 1, 1, "ZZZZ", "0");
}

#[test]
fn sim1_3() {
    run_test("sim1", 4, 2, 5, "ZZZZ", "2^-12 * (27)");
}

#[test]
fn sim2_1() {
    run_test("sim2", 6, 2, 0, "ZZZZZZ", "2^-14 * (243)");
}

#[test]
fn sim2_2() {
    run_test("sim2", 4, 2, 9, "ZXXZ", "2^-12 * (273)");
}

#[test]
fn sim9_1() {
    run_test("sim9", 5, 3, 0, "ZZZZZ", "2^-13 * (131)");
}
