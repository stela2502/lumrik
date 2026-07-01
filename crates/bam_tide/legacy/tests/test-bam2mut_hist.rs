use std::process::Command;
use std::process::exit;

#[test]
fn test_bam2mut_hist_output() {

    let is_release_mode = !cfg!(debug_assertions);


    let exe = env!("CARGO_BIN_EXE_bam-coverage");

    print!("{}",command);

    let output = Command::new(command)
        .args(&["--input", "testData/bowtie_chrm_mutated.bam", "-b", "10", "-H"])
        .output()
        .expect("Failed to run bam2mut_hist");

    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);

    let expected: Vec<usize> = vec![1904, 1416, 1422, 1504, 1402, 1452, 1399, 1395, 1671, 1603];

    // Collect every parsed "index : value" pair into a Vec<[i32; 2]>
    //print!("'{stdout}'");
    let parsed: Vec<[usize; 2]> = stdout
        .lines()
        .filter_map(|line| {
            let mut parts = line.splitn(2, ':');
            let idx_str = parts.next()?.trim();
            let num_str = parts.next()?.trim();

            // Parse both pieces as i32 and pack into a fixed-size array
            match (idx_str.parse::<usize>(), num_str.parse::<usize>()) {
                (Ok(idx), Ok(num)) => Some([idx, num]),
                _ => None,
            }
        })
        .collect();
    for (i, &val) in expected.iter().enumerate() {
        // Find the line that starts with the index
        
        assert_eq!(parsed[i][0], i, "Wrong index in line: {i}");
        assert_eq!(parsed[i][1], val, "Wrong value in line: {i}");
    }
}