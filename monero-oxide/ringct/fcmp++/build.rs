use std::{
  io::Write,
  env,
  path::Path,
  fs::{File, remove_file},
};

#[cfg(feature = "compile-time-generators")]
trait NewGenerators: ciphersuite::Ciphersuite {
  fn new_generators() -> monero_generators::FcmpGenerators<Self>;
}
#[cfg(feature = "compile-time-generators")]
impl NewGenerators for helioselene::Helios {
  fn new_generators() -> monero_generators::FcmpGenerators<Self> {
    monero_generators::FcmpGenerators::<Self>::new()
  }
}
#[cfg(feature = "compile-time-generators")]
impl NewGenerators for helioselene::Selene {
  fn new_generators() -> monero_generators::FcmpGenerators<Self> {
    monero_generators::FcmpGenerators::<Self>::new()
  }
}

#[cfg(feature = "compile-time-generators")]
fn generator_set<C: NewGenerators>()
where
  C::G: ciphersuite::group::GroupEncoding<Repr = [u8; 32]> + ec_divisors::DivisorCurve,
{
  use ciphersuite::group::{ff::PrimeField, GroupEncoding};
  use ec_divisors::DivisorCurve;

  fn serialize<G: GroupEncoding<Repr = [u8; 32]> + DivisorCurve>(
    generators_string: &mut String,
    point: &str,
    points: &[G],
  ) {
    for generator in points {
      let (x, y) = G::to_xy(*generator).expect("generator was the identity?");
      generators_string.extend(
        format!(
          "
            helioselene::{point}::from_xy(
              <
                <helioselene::{point} as ec_divisors::DivisorCurve>::FieldElement
                  as
                ciphersuite::group::ff::PrimeField
              >::from_repr({:?}).expect(\"build script x wasn't reduced\"),
              <
                <helioselene::{point} as ec_divisors::DivisorCurve>::FieldElement
                  as
                ciphersuite::group::ff::PrimeField
              >::from_repr({:?}).expect(\"build script y wasn't reduced\"),
            ).expect(\"generator from build script wasn't on-curve\"),
          ",
          x.to_repr().as_ref(),
          y.to_repr().as_ref()
        )
        .chars(),
      );
    }
  }

  let id =
    String::from_utf8(C::ID.to_vec()).expect("Helios/Selene Ciphersuite ID wasn't valid UTF-8");
  let point = format!("{id}Point");
  let path = format!("{}_generators.rs", id.clone().to_lowercase());

  let generators = C::new_generators();
  let mut g = String::new();
  serialize(&mut g, &point, &[generators.generators.g()]);
  let mut h = String::new();
  serialize(&mut h, &point, &[generators.generators.h()]);

  let mut g_bold: Vec<u8> =
    Vec::with_capacity((2 * 32) * generators.generators.g_bold_slice().len());
  for generator in generators.generators.g_bold_slice() {
    let (x, y) = C::G::to_xy(*generator).expect("generator was the identity?");
    g_bold.extend(x.to_repr().as_ref());
    g_bold.extend(y.to_repr().as_ref());
  }

  let mut h_bold: Vec<u8> =
    Vec::with_capacity((2 * 32) * generators.generators.h_bold_slice().len());
  for generator in generators.generators.h_bold_slice() {
    let (x, y) = C::G::to_xy(*generator).expect("generator was the identity?");
    h_bold.extend(x.to_repr().as_ref());
    h_bold.extend(y.to_repr().as_ref());
  }

  let path = Path::new(&env::var("OUT_DIR").expect("cargo didn't set $OUT_DIR")).join(path);
  let _ = remove_file(&path);
  File::create(&path)
    .expect("failed to create file in $OUT_DIR")
    .write_all(
      format!(
        "
          /// The FCMP generators for {id}.
          pub(super) static {}_FCMP_GENERATORS:
            std_shims::sync::LazyLock<monero_generators::FcmpGenerators<helioselene::{}>> =
              std_shims::sync::LazyLock::new(|| monero_generators::FcmpGenerators {{
                generators: generalized_bulletproofs::Generators::new(
                  {g}
                  {h}
                  {{
                    const BYTES: &[u8] = &{:?};
                    let mut bytes = BYTES;
                    let mut x = [0; 32];
                    let mut y = [0; 32];
                    let mut res = Vec::with_capacity({});
                    while !bytes.is_empty() {{
                      x.copy_from_slice(&bytes[.. 32]);
                      bytes = &bytes[32 ..];
                      y.copy_from_slice(&bytes[.. 32]);
                      bytes = &bytes[32 ..];
                      res.push(
                        helioselene::{point}::from_xy(
                          <
                            <helioselene::{point} as ec_divisors::DivisorCurve>::FieldElement
                              as
                            ciphersuite::group::ff::PrimeField
                          >::from_repr(x).expect(\"build script x wasn't reduced\"),
                          <
                            <helioselene::{point} as ec_divisors::DivisorCurve>::FieldElement
                              as
                            ciphersuite::group::ff::PrimeField
                          >::from_repr(y).expect(\"build script y wasn't reduced\"),
                      ).expect(\"generator from buildf script wasn't on-curve\")
                      );
                    }}
                    res
                  }},
                  {{
                    const BYTES: &[u8] = &{:?};
                    let mut bytes = BYTES;
                    let mut x = [0; 32];
                    let mut y = [0; 32];
                    let mut res = Vec::with_capacity({});
                    while !bytes.is_empty() {{
                      x.copy_from_slice(&bytes[.. 32]);
                      bytes = &bytes[32 ..];
                      y.copy_from_slice(&bytes[.. 32]);
                      bytes = &bytes[32 ..];
                      res.push(
                        helioselene::{point}::from_xy(
                          <
                            <helioselene::{point} as ec_divisors::DivisorCurve>::FieldElement
                              as
                            ciphersuite::group::ff::PrimeField
                          >::from_repr(x).expect(\"build script x wasn't reduced\"),
                          <
                            <helioselene::{point} as ec_divisors::DivisorCurve>::FieldElement
                              as
                            ciphersuite::group::ff::PrimeField
                          >::from_repr(y).expect(\"build script y wasn't reduced\"),
                        ).expect(\"generator from build script wasn't on-curve\")
                      );
                    }}
                    res
                  }},
                ).unwrap()
              }});
        ",
        id.clone().to_uppercase(),
        id,
        g_bold.as_slice(),
        generators.generators.g_bold_slice().len(),
        h_bold.as_slice(),
        generators.generators.h_bold_slice().len(),
      )
      .as_bytes(),
    )
    .expect("couldn't write generated source code to file on disk");
}

#[cfg(feature = "compile-time-generators")]
fn generators() {
  generator_set::<helioselene::Helios>();
  generator_set::<helioselene::Selene>();
  let path =
    Path::new(&env::var("OUT_DIR").expect("cargo didn't set $OUT_DIR")).join("generators.rs");
  let _ = remove_file(&path);
  File::create(&path)
    .expect("failed to create file in $OUT_DIR")
    .write_all(
      b"
        mod helios_generators {
          include!(concat!(env!(\"OUT_DIR\"), \"/helios_generators.rs\"));
        }
        mod selene_generators {
          include!(concat!(env!(\"OUT_DIR\"), \"/selene_generators.rs\"));
        }
        use helios_generators::HELIOS_FCMP_GENERATORS;
        use selene_generators::SELENE_FCMP_GENERATORS;
      ",
    )
    .expect("couldn't write generated source code to file on disk");
}

#[cfg(not(feature = "compile-time-generators"))]
fn generators() {
  let path =
    Path::new(&env::var("OUT_DIR").expect("cargo didn't set $OUT_DIR")).join("generators.rs");
  let _ = remove_file(&path);
  File::create(&path)
    .expect("failed to create file in $OUT_DIR")
    .write_all(
      format!(
        "
          /// The FCMP generators for Helios.
          static HELIOS_FCMP_GENERATORS: LazyLock<monero_generators::FcmpGenerators<Helios>> =
            LazyLock::new(monero_generators::FcmpGenerators::<Helios>::new);
          /// The FCMP generators for Selene.
          static SELENE_FCMP_GENERATORS: LazyLock<monero_generators::FcmpGenerators<Selene>> =
            LazyLock::new(monero_generators::FcmpGenerators::<Selene>::new);
      ",
      )
      .as_bytes(),
    )
    .expect("couldn't write generated source code to file on disk");
}

fn main() {
  println!("cargo:rerun-if-changed=build.rs");

  generators();
}
