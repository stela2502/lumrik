/// Removes `option` and its following arguments.
///
/// `Some(n)` removes exactly `n` arguments after the option.
/// `None` removes arguments until the next `-` or `--` option.
pub fn remove_option(args: &mut Vec<String>, option: &str, n: Option<usize>) {
    while let Some(pos) = args.iter().position(|arg| arg == option) {
        let end = match n {
            Some(n) => (pos + 1 + n).min(args.len()),
            None => args[pos + 1..]
                .iter()
                .position(|arg| arg.starts_with('-'))
                .map(|next| pos + 1 + next)
                .unwrap_or(args.len()),
        };

        args.drain(pos..end);
    }
}

pub fn has_option(args: &[String], option: &str) -> bool {
    args.iter().any(|arg| arg == option)
}
