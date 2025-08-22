use std::{
  io::Write,
  env,
  path::Path,
  fs::{File, remove_file},
};

#[cfg(feature = "compile-time-generators")]
fn generators() {
  use ciphersuite::group::GroupEncoding;
  use helioselene::{Helios, Selene};
  use monero_generators::FcmpGenerators;

  fn serialize<G: GroupEncoding<Repr = [u8; 32]>>(
    generators_string: &mut String,
    point: &'static str,
    points: &[G],
  ) {
    for generator in points {
      generators_string.extend(
        format!(
          "
            helioselene::{point}::from_bytes({:?})
              .expect(\"generator from build script wasn't on-curve\"),
          ",
          generator.to_bytes()
        )
        .chars(),
      );
    }
  }

  let helios = FcmpGenerators::<Helios>::new();
  let mut helios_g = String::new();
  serialize(&mut helios_g, "HeliosPoint", &[helios.generators.g()]);
  let mut helios_h = String::new();
  serialize(&mut helios_h, "HeliosPoint", &[helios.generators.h()]);
  let mut helios_g_bold = String::new();
  serialize(&mut helios_g_bold, "HeliosPoint", helios.generators.g_bold_slice());
  let mut helios_h_bold = String::new();
  serialize(&mut helios_h_bold, "HeliosPoint", helios.generators.h_bold_slice());

  let selene = FcmpGenerators::<Selene>::new();
  let mut selene_g = String::new();
  serialize(&mut selene_g, "SelenePoint", &[selene.generators.g()]);
  let mut selene_h = String::new();
  serialize(&mut selene_h, "SelenePoint", &[selene.generators.h()]);
  let mut selene_g_bold = String::new();
  serialize(&mut selene_g_bold, "SelenePoint", selene.generators.g_bold_slice());
  let mut selene_h_bold = String::new();
  serialize(&mut selene_h_bold, "SelenePoint", selene.generators.h_bold_slice());

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
            LazyLock::new(monero_generators::FcmpGenerators {{
              generators: generalized_bulletproofs::Generators::new(
                {helios_g}
                {helios_h}
                vec![{helios_g_bold}],
                vec![{helios_h_bold}],
              )
            }});

          /// The FCMP generators for Selene.
          static SELENE_FCMP_GENERATORS: LazyLock<monero_generators::FcmpGenerators<Selene>> =
            LazyLock::new(monero_generators::FcmpGenerators {{
              generators: generalized_bulletproofs::Generators::new(
                {selene_g}
                {selene_h}
                vec![{selene_g_bold}],
                vec![{selene_h_bold}],
              )
            }});
        ",
      )
      .as_bytes(),
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
