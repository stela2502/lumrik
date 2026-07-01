#! bash

find target/x86_64-unknown-linux-musl/release/ \
	-maxdepth 1 \
	-type f \
	-executable \
	-exec cp {} /home/med-sal/sens05_home/bin/ \;
