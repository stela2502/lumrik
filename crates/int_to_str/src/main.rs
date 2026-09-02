use clap::Parser;
use int_to_str::int_to_str::IntToStr;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Convert DNA sequences to Lumrik integer encoding and back"
)]
struct Cli {
    /// DNA sequence or encoded unsigned integer.
    #[arg(value_name = "SEQUENCE_OR_INTEGER")]
    input: String,
}

fn main() {
    let cli = Cli::parse();

    if let Ok(num) = cli.input.parse::<u128>() {
        let tool = IntToStr::from_u128(num);
        println!("Integer input: {num}");
        println!("→ Sequence: {}", tool.to_string(64));
    } else {
        println!("Sequence input: {}", cli.input);
        let tool = IntToStr::new(&cli.input);
        println!("→ Sequence: {}", tool.into_u128());
    }
}
