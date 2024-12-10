use circuit_runner::converter::convert_poc_graph_to_quizx;
use quizx::decompose_tri::*;
use quizx::graph::*;
use quizx::hash_graph::*;
use std::env;
use std::path::Path;
use std::time::Instant;

fn main() {
    let args: Vec<String> = env::args().collect();

    let without_arg = args.get(1).expect("Missing 'without' argument");
    let with_sd = match without_arg.as_str() {
        "WITH_SD" => true,
        "WITHOUT_SD" => false,
        _ => panic!("Invalid value for 'without' argument, expected 'WITH_SD' or 'WITHOUT_SD'"),
    };

    let file_path = args.get(2).expect("Missing file path");
    let file_path = Path::new(file_path);

    let dot_string = std::fs::read_to_string(file_path).expect("Failed to read file");
    let dot_string = dot_string.replace("\r", "");
    // println!("received string:\n{}", dot_string);

    let mut poc_graph: Graph = GraphLike::from_dot(&dot_string);

    // println!("Graph loaded!");

    let mut quizx_graph = convert_poc_graph_to_quizx(&mut poc_graph);

    // println!(
    //     "----------------------------------------------------------------------------------------"
    // );
    // graphviz-friendly string representation of the graph
    // println!("string of final diagram:\n{}", quizx_graph.to_dot());

    let time = Instant::now();

    quizx::simplify::clifford_simp(&mut quizx_graph);

    // println!("{}", time.elapsed().as_secs_f64());

    let mut d = TriDecomposer::new();
    d.add_decomp(Decompositions::ZSpider(ZSpiderDecomp::new()));
    if with_sd {
        d.add_decomp(Decompositions::Tri(TriDecomp::new(3)));
    }

    let (_scalar, num_terms) = d.decomp_sequential(quizx_graph);

    // println!("Scalar: {}", scalar);

    println!("QUIZX_WITH NUM OF TERMS: {}", num_terms); // oopsie this naming should have been dependent on with_sd variable

    println!("QUIZX_WITH RUNTIME: {}", time.elapsed().as_secs_f64());
}
