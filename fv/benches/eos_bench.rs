extern crate fv;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use fv::eos::eos::EquationOfState;

fn bench_conserved_to_primitive(b: &mut Criterion) {
    let mut eos = EquationOfState::new(1000, 1000);
    eos.u.data.fill(1.0);

    b.bench_function("conserved_to_primitive", |b| {
        b.iter(|| eos.conserved_to_primitive())
    });
}

criterion_group!(benches, bench_conserved_to_primitive);
criterion_main!(benches);
