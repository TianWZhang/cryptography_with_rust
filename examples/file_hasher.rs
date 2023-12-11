use cryptography_with_rust::sha256::{self, Sha256};
use std::{
    env,
    fs::File,
    io::{BufRead, BufReader},
};

fn main() -> std::io::Result<()> {
    let mut hasher = Sha256::new();
    let args = env::args().collect::<Vec<String>>();

    if args.len() != 2 {
        println!("Usage:");
        println!("file_hasher: <filename>");
        return Ok(());
    }

    let filename = args[1].trim();
    let file = File::open(filename)?;
    let chunk_size = sha256::BLOCK_SIZE * 1024;
    let mut reader = BufReader::with_capacity(chunk_size, &file);

    // read a file in byte chunks
    loop {
        let buffer = reader.fill_buf()?;
        hasher.update(buffer);
        let len = buffer.len();
        if len == 0 {
            break;
        }
        reader.consume(len);
    }

    println!("{} {}", hasher.hex_digest().to_ascii_lowercase(), filename);
    Ok(())
}
