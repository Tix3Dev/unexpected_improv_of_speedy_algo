use std::env;
use std::time::Instant;

use quizx::circuit::Circuit;
use quizx::decompose_tri::*;
use zxbarren::variance::*;
use zxbarren::sim;
use zxbarren::tree;

fn get_ansatz(name: &String, n: usize, l: usize) -> (Circuit, Vec<String>) {
    match name.as_str() {
        "sim1" => sim::circ1(n, l),
        "sim2" => sim::circ2(n, l),
        "sim9" => sim::circ9(n, l),
        "sim10" => sim::circ10(n, l),
        "sim11" => sim::circ11(n, l),
        "sim12" => sim::circ12(n, l),
        "sim15" => sim::circ15(n, l),
        "tree" => tree::tree(n),
        _ => panic!("Circuit not found: {name}")
    }
}



fn main() {
    const QUIZX: &str = "quizx";
    const OURS: &str = "ours";

    let args: Vec<_> = env::args().collect();
    let (name, n, l, mode, parallel) = (args[1].parse().unwrap(),
                                args[2].parse().unwrap(),
                                args[3].parse().unwrap(),
                                args[4].parse::<String>().unwrap(),
                                args[5].parse().unwrap());

        
    let (c, params) = get_ansatz(&name, n, l);

    let time = Instant::now();

    if mode == QUIZX {
        let mut var = create_diagram(&c, &"Z".repeat(n), &params[0], false);
        quizx::simplify::full_simp(&mut var);
        let mut d = quizx::decompose::Decomposer::new(&var);
        d.with_full_simp();
        d.use_cats(true);
        if parallel {
            d.decomp_parallel(3);
        } else {
            d.decomp_all();
        }
    } else if mode == OURS {
        let var_orig = create_diagram(&c, &"Z".repeat(n), &params[0], true);

        let mut var = var_orig.clone();
        quizx::simplify::clifford_simp(&mut var);
        let mut d = TriDecomposer::new();
        d.add_decomp(Decompositions::ZSpider(ZSpiderDecomp::new()));
        d.add_decomp(Decompositions::Tri(TriDecomp::new(3)));
        // if parallel {
        //     d.decomp_parallel(var);
        // } else {
        //     d.decomp_sequential(var);
        // }
        let (_scalar, num_terms) = d.decomp_sequential(var);
        let time1 = time.elapsed().as_secs_f64();
        println!("\n{}", time1);
        println!("QUIZX_WITH NUM OF TERMS: {}", num_terms);


        let mut var = var_orig.clone();
        quizx::simplify::clifford_simp(&mut var);
        let mut d = TriDecomposer::new();
        d.add_decomp(Decompositions::ZSpider(ZSpiderDecomp::new()));
        // d.add_decomp(Decompositions::Tri(TriDecomp::new(3)));
        // if parallel {
        //     d.decomp_parallel(var);
        // } else {
        //     d.decomp_sequential(var);
        // }
        let (_scalar, num_terms) = d.decomp_sequential(var);
        let time2 = time.elapsed().as_secs_f64() - time1;
        println!("{}", time2);
        println!("QUIZX_WITHOUT NUM OF TERMS: {}\n", num_terms);
    }
    if mode == QUIZX { // QUIZX is called after OURS, so this will end one comparison
        println!("\n{}", time.elapsed().as_secs_f64());
        println!("---------------------------------------------------------------------------");
    }
}
