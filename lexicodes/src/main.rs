use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use std::collections::HashSet;

fn binomial(n: usize, k: usize) -> usize {
    (0..k).fold(1, |acc, i| acc * (n - i) / (i + 1))
}

fn to_bitmask(comb: &[usize]) -> u32 {
    comb.iter().fold(0, |acc, &x| acc | (1 << x))
}

fn lexicode(k: usize) -> HashSet<u32> {
    let mut added: HashSet<u32> = HashSet::new();
    let total_combinations = binomial(2 * k, k - 1) as u64;

    let pb = ProgressBar::new(total_combinations);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );

    let chunk_size = 10_000; // Update the bar every 10,000 iterations
    for (i, c) in (0..2 * k).combinations(k - 1).enumerate() {
        let c_mask = to_bitmask(&c);
        let mut valid = true;
        for i in 0..2 * k {
            for j in 0..i {
                if ((c_mask >> i) ^ (c_mask >> j)) & 1 == 1 {
                    let altered = c_mask ^ (1 << i) ^ (1 << j);
                    if added.contains(&altered) {
                        valid = false;
                        break;
                    }
                }
            }
            if !valid {
                break;
            }
        }
        if valid {
            added.insert(c_mask);
        }
        if i % chunk_size == 0 {
            pb.set_position(i as u64);
        }
    }
    pb.finish_and_clear();

    added
}

fn check_intersections(k: u32, added: &HashSet<u32>) -> bool {
    let pb = ProgressBar::new(added.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );

    let largest: u32 = (1u32.checked_shl(2 * k).unwrap_or(0)).saturating_sub(1);
    for ai in added.iter() {
        for i in 0..2 * k {
            let ai_complement = ai ^ largest;
            if ai_complement >> i & 1 != 1 {
                continue;
            }
            for j in 0..i {
                if ai_complement >> j & 1 == 1 {
                    let ij_removed = ai_complement ^ (1 << i) ^ (1 << j);
                    if added.contains(&ij_removed) {
                        pb.finish_and_clear();
                        return false;
                    }
                }
            }
        }
        pb.inc(1);
    }
    pb.finish_and_clear();
    return true;
}

fn main() {
    for k in 2..=16 {
        println!("Building lexicode for k = {}", k);
        let added = lexicode(k);
        let num_words = added.len();
        let target = binomial(2 * k, k) / (2 * (k + 1));
        println!("Found {} codewords. Target: {}.", num_words, target);
        if num_words == target {
            println!("Target match; checking intersections.");
            let is_good = check_intersections(k as u32, &added);
            if is_good {
                println!(
                    "Intersections good, θ(𝐽({}, {})) = {}",
                    2 * k,
                    k,
                    binomial(2 * k, k) / (k + 1)
                );
            } else {
                println!("Found empty intersection.");
            }
        }
        println!("");
    }
}
