#![no_std]

#[cfg(feature = "alloc")]
pub mod alloc {
  pub use helioselene;
  pub use monero_wallet;
}
