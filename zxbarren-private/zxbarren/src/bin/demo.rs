use zxbarren::sim;
use zxbarren::variance::*;

fn main() {
    for n in 1..33 {
        let h = "Z".repeat(n);
        let (c, params) = sim::circ10(n, 1);

        let var = compute_variance(&c, &h, &params[0]);
        println!("{n} : {}", var.to_float());
    }
}
