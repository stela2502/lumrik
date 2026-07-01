use std::fs;
//use std::collections::HashMap;
use std::path::Path;
use std::io::BufRead;
use std::fs::File;
use std::io;

use bam_tide::bed_data::BedData;
use bam_tide::gtf_logics::{AnalysisType};

#[test]

fn test_bam_to_bigwig() {
    let bam_file = "testData/test.bam"; // Replace with a small BAM test file
    let bigwig_file = "testData/output/test.wig";

    // Ensure output folder exists
    fs::create_dir_all("testData/output/").unwrap();

    let path = Path::new(bigwig_file);
    if path.exists() {
        if let Err(e) = fs::remove_file(path) {
            eprintln!("Failed to remove file {}: {}", bigwig_file, e);
        }
    }

    // Run the function
    let data =  BedData::new( bam_file, 50, 1, &AnalysisType::Bulk, &b"CR", &b"No", false, false, 1_u8 );

    let result = BedData::write_bedgraph( &data, bigwig_file );

    // Assert success
    assert!(result.is_ok());

    // Optionally, validate the BigWig file
    assert!(fs::metadata(bigwig_file).is_ok());
}



#[test]
fn test_bam_to_bedgraph() {
    let bam_file = "testData/minimal_test.bam"; // Replace with a small BAM test file
    let bg_file = "testData/output/mini_test.bg";

    // Ensure output folder exists
    fs::create_dir_all("testData/output/").unwrap();

    let path = Path::new(bg_file);
    if path.exists() {
        if let Err(e) = fs::remove_file(path) {
            eprintln!("Failed to remove file {}: {}", bg_file, e);
        }
    }

    // Run the function
    let data =  BedData::new( bam_file, 50, 1, &AnalysisType::Bulk, &b"CR", &b"No", false, false, 1_u8 );
    let result = BedData::write_bedgraph( &data, bg_file );

    // Assert success
    assert!(result.is_ok());

    // Optionally, validate the BigWig file
    assert!(fs::metadata(bg_file).is_ok());

    

    let exp = vec![
        ["chr1".to_string(), "0".to_string(), "14400".to_string(), "0".to_string()],
        ["chr1".to_string(), "14400".to_string(), "14450".to_string(), "1".to_string()],
        ["chr1".to_string(), "14450".to_string(), "14500".to_string(), "3".to_string()],
        ["chr1".to_string(), "14500".to_string(), "14550".to_string(), "8".to_string()],
        ["chr1".to_string(), "14550".to_string(), "14600".to_string(), "12".to_string()],
        ["chr1".to_string(), "14600".to_string(), "14650".to_string(), "10".to_string()],
        ["chr1".to_string(), "14650".to_string(), "14700".to_string(), "4".to_string()],
        ["chr1".to_string(), "14700".to_string(), "14850".to_string(), "1".to_string()],
        ["chr1".to_string(), "14850".to_string(), "16600".to_string(), "0".to_string()],
    ];
    match read_tsv_to_vec( bg_file ) {
        Ok(data) => {
            let first_10 = &data[..10.min(data.len())];
            for (i, (expected, actual)) in exp.iter().zip( first_10.iter()).enumerate() {
                assert_eq!(actual, expected, "Mismatch at line {}", i + 1);
            }
        },
        Err(e) => panic!("We reached an unexpected error reading the ofile {bg_file}: {e:?}"),
    }
}

fn read_tsv_to_vec(path: &str) -> io::Result<Vec<[String; 4]>> {
    let file = File::open(path)?;
    let reader = io::BufReader::new(file);
    
    let mut result = Vec::new();
    
    for line in reader.lines() {
        let line = line?; // Handle potential I/O error
        let fields: Vec<String> = line.split('\t').map(String::from).collect();
        
        if fields.len() == 4 {
            result.push([fields[0].clone(), fields[1].clone(), fields[2].clone(), fields[3].clone()]);
        } else {
            eprintln!("Skipping invalid line: {}", line);
        }
    }
    
    Ok(result)
}