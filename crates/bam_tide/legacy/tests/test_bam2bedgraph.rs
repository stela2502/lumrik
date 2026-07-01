// Import necessary modules
use std::process::Command;
//use std::fs::File;
//use std::fs;
//use std::io::BufReader;
//use std::io::BufRead;
//use std::collections::HashMap;
use std::process::exit;

// bigwig should work as well - but it creates enormouse files as it adds the zeros in, too.
#[test]
fn test_bam2bedgraph() {

    let is_release_mode = !cfg!(debug_assertions);

    if ! is_release_mode {
        eprintln!("Test should be re-run in release mode (speed!)");
        return;
    }

    let command = env!("CARGO_BIN_EXE_bam2bedgraph");

    let args = &[
        "-b", "testData/test.bam",
        "-o", "testData/output/test.bedgraph",
        "-w", "50", 
        "-u", "No",
        "-a", "bulk",
    ]; 


    // Execute the command with the provided arguments
    let output = Command::new( command ).args( args )
        .output()
        .map_err(|e| {
            eprintln!("Failed to execute command: {}", e);
            e
        }).unwrap();

    let cmd = format!("{} {}", command, args.join(" "));
    if !output.status.success() {
        eprintln!("Command failed: {}", cmd);
        // Handle failure accordingly
    }else {
        println!("{}", cmd );
    }

    // Check if the command was successful (exit code 0)
    assert!(output.status.success());
    // Convert output to string
    let _output_str = String::from_utf8_lossy(&output.stdout);

}