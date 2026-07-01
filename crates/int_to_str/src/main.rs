use std::env;
use stela_int_to_str::int_to_str::IntToStr;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        eprintln!("Usage: {} <sequence or integer>", args[0]);
        return;
    }

    let input = &args[1];

    // Try parsing as an integer
    if let Ok(num) = input.parse::<u128>() {
        let tool = IntToStr::from_u128(num);

        println!("Integer input: {}", num);
        let seq = tool.to_string(64);
        println!("→ Sequence: {}", seq);
    } else {
        println!("Sequence input: {}", input);
        let tool = IntToStr::new(input);
        println!("→ Sequence: {}", tool.into_u128());
    }
}
