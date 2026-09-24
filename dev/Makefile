# Makefile for backing up a GitHub repository using rsync

SOURCE_DIR := ./Rustody
DEST_DIR := stefanl@aurora-ls2.lunarc.lu.se:~/NAS/git/Rustody

backup:
	rsync -av --exclude='.git*' $(SOURCE_DIR)/ $(DEST_DIR)/

Rustody_to_LS2: backup
clean:
	@for path in */Cargo.toml; do \
		if [ -f "$$path" ]; then \
			dir=$$(dirname "$$path"); \
			echo "Cleaning in $$dir"; \
			(cd "$$dir" && cargo clean); \
		fi \
	done

RUST_WORKFLOW := .github/workflows/rust.yml

rust-setup:
	@echo "Scanning for Cargo.toml files..."
	@find ./ -maxdepth 2 -name Cargo.toml | while read cargo; do \
		dir=$$(dirname $$cargo); \
		repo_name=$$(basename $$dir); \
		readme="$$dir/README.md"; \
		workflow="$$dir/$(RUST_WORKFLOW)"; \
		\
		rust_badge="[![Rust](https://github.com/stela2502/$$repo_name/actions/workflows/rust.yml/badge.svg)](https://github.com/stela2502/$$repo_name/actions/workflows/rust.yml)"; \
		\
		if [ ! -f "$$readme" ]; then \
			echo "Creating $$readme"; \
			echo "$$rust_badge" > "$$readme"; \
		else \
			first_line=$$(head -n1 "$$readme"); \
			if [ "$$first_line" != "$$rust_badge" ]; then \
				echo "Fixing $$readme first line"; \
				tmp=$$(mktemp); \
				echo "$$rust_badge" > $$tmp; \
				cat "$$readme" >> $$tmp; \
				mv $$tmp "$$readme"; \
			fi \
		fi; \
		\
		if [ ! -f "$$workflow" ]; then \
			echo "Copying workflow to $$workflow"; \
			mkdir -p "$$(dirname $$workflow)"; \
			cp rust_cli.yml "$$workflow"; \
		fi; \
	done

.PHONY: install-musl-bins lumrik_tar_file

ROOT := $(HOME)/git_Projects
MUSL_TARGET := target/x86_64-unknown-linux-musl/release

DEST1 := $(HOME)/sens05_home/bin
DEST2 := $(HOME)/sens05_shared/common/software/Rustody/1.5/bin

lumrik_tar_file:
	tar -czf lumrik-source.tar.gz \
		--exclude='*/target' \
		--exclude='*/target/*' \
		--exclude='*/out/*' \
		--exclude='*/.git' \
		--exclude='*/.git/*' \
		--exclude='*.tar.gz' \
		--exclude='_STARtmp' \
		--exclude='_STARtmp/*' \
		--exclude='nelrune_test_out' \
		--exclude='nelrune_test_out/*' \
		--exclude='*.bam' \
		--exclude='*.bai' \
		--exclude='*.cram' \
		--exclude='*.fastq' \
		--exclude='*.fastq.gz' \
		--exclude='*.fq' \
		--exclude='*.fq.gz' \
		--exclude='*.sif' \
		--exclude='*.rds' \
		--exclude='*.h5' \
		--exclude='*.h5ad' \
		--exclude='*.mtx' \
		--exclude='*.mtx.gz' \
		--exclude='*/whitelists/*.txt.gz' \
		--exclude='*/star_index' \
		--exclude='*/.idea' \
		--exclude='*/.vscode' \
		--exclude='*.swp' \
		--exclude='sequences_many.txt' \
		lumrik
	@du -h lumrik-source.tar.gz


MUSL_TRIPLE := x86_64-unknown-linux-musl
MUSL_TARGET := $(LUMRIK_DIR)/target/$(MUSL_TRIPLE)/release

install-musl-bins:
	@mkdir -p "$(DEST1)" "$(DEST2)" "$(HOME)/bin"
	@echo "Building lumrik musl release binaries"
	@cd "$(LUMRIK_DIR)" && cargo build -r --target $(MUSL_TRIPLE)
	@echo "Scanning $(LUMRIK_DIR)/$(MUSL_TARGET)"
	@find "$(LUMRIK_DIR)/$(MUSL_TARGET)" \
		-maxdepth 1 \
		-type f \
		-executable \
		! -name "*.d" \
		-print | while read bin; do \
			name=$$(basename "$$bin"); \
			for dest in "$(DEST1)" "$(DEST2)" "$(HOME)/bin"; do \
				target="$$dest/$$name"; \
				if [ ! -e "$$target" ] || ! cmp -s "$$bin" "$$target"; then \
					echo "  installing $$name -> $$dest"; \
					cp "$$bin" "$$target"; \
					chmod 755 "$$target"; \
				else \
					echo "  unchanged  $$name -> $$dest"; \
				fi; \
			done; \
		done

install-local:
	@echo "Installing musl binaries into /usr/local/bin"
	@find "$(ROOT)" -maxdepth 2 -name Cargo.toml -print | while read cargo; do \
		crate_dir=$$(dirname "$$cargo"); \
		release_dir="$$crate_dir/$(MUSL_TARGET)"; \
		if [ ! -d "$$release_dir" ]; then \
			continue; \
		fi; \
		echo "Scanning $$crate_dir"; \
		find "$$release_dir" -maxdepth 1 -type f -executable ! -name "*.d" -print | while read bin; do \
			name=$$(basename "$$bin"); \
			target="/usr/local/bin/$$name"; \
			if [ ! -e "$$target" ] || ! cmp -s "$$bin" "$$target"; then \
				echo "  installing $$name -> /usr/local/bin"; \
				sudo install -m 755 "$$bin" "$$target"; \
			else \
				echo "  unchanged  $$name -> /usr/local/bin"; \
			fi; \
		done; \
	done

install-sc-ont-mut-quant:
	make -C sc-ont-mut-quant clean
	rsync -av \
		--no-owner \
		--no-group \
		--exclude '.git' \
		--exclude '.github' \
		--exclude '.build' \
		--exclude '.bam' \
		--exclude 'work' \
		--exclude '.nextflow*' \
		sc-ont-mut-quant/ \
		~/sens05_shared/nf-core/sc-ont-mut-quant/
LUMRIK_DIR := $(ROOT)/lumrik
LUMRIK_IMAGE := ghcr.io/stela2502/lumrik
SINGULARITY_LOCAL := /data1/nextflow/singularity
SINGULARITY_SENS := $(HOME)/sens05_home/nextflow/singularity_cache

LUMRIK_DIR := $(ROOT)/lumrik
LUMRIK_IMAGE := ghcr.io/stela2502/lumrik
SINGULARITY_LOCAL := /data1/nextflow/singularity
SINGULARITY_SENS := $(HOME)/sens05_home/nextflow/singularity_cache

lumrik-release-image:
	@set -e; \
	cd "$(LUMRIK_DIR)"; \
	if [ -n "$$(git status --porcelain)" ]; then \
		echo "ERROR: lumrik repository contains uncommitted changes"; \
		git status --short; \
		exit 1; \
	fi; \
	git fetch --tags --quiet; \
	last_tag=$$(git tag --list 'v[0-9]*.[0-9]*.[0-9]*' --sort=-v:refname | head -n1); \
	requested_version="$(VERSION)"; \
	if [ -n "$$requested_version" ]; then \
		version="$$requested_version"; \
	else \
		if [ -z "$$last_tag" ]; then \
			version="0.1.0"; \
		else \
			previous=$${last_tag#v}; \
			major=$${previous%%.*}; \
			rest=$${previous#*.}; \
			minor=$${rest%%.*}; \
			patch=$${rest##*.}; \
			patch=$$((patch + 1)); \
			version="$$major.$$minor.$$patch"; \
		fi; \
	fi; \
	if ! printf '%s\n' "$$version" | grep -Eq '^[0-9]+\.[0-9]+\.[0-9]+$$'; then \
		echo "ERROR: VERSION must be X.Y.Z, got: $$version"; \
		exit 1; \
	fi; \
	tag="v$$version"; \
	echo "Latest tag:      $${last_tag:-none}"; \
	echo "Release version: $$version"; \
	echo "Release tag:     $$tag"; \
	if git rev-parse -q --verify "refs/tags/$$tag" >/dev/null; then \
		tag_commit=$$(git rev-list -n1 "$$tag"); \
		head_commit=$$(git rev-parse HEAD); \
		if [ "$$tag_commit" != "$$head_commit" ]; then \
			echo "ERROR: $$tag already exists and does not point at HEAD"; \
			exit 1; \
		fi; \
		cargo_version=$$(awk '/^\[workspace.package\]$$/{p=1; next} /^\[/{p=0} p && /^version = "/{gsub(/^version = "|"$$/, ""); print; exit}' Cargo.toml); \
		if [ "$$cargo_version" != "$$version" ]; then \
			echo "ERROR: $$tag exists at HEAD but Cargo.toml says $$cargo_version"; \
			exit 1; \
		fi; \
		echo "Existing release $$tag found at HEAD; reusing it."; \
	else \
		old_version=$$(awk '/^\[workspace.package\]$$/{p=1; next} /^\[/{p=0} p && /^version = "/{gsub(/^version = "|"$$/, ""); print; exit}' Cargo.toml); \
		if [ -z "$$old_version" ]; then \
			echo "ERROR: could not read [workspace.package] version from Cargo.toml"; \
			exit 1; \
		fi; \
		if [ "$$old_version" = "$$version" ]; then \
			echo "ERROR: Cargo.toml is already $$version but tag $$tag does not exist"; \
			exit 1; \
		fi; \
		echo "Updating Cargo workspace version: $$old_version -> $$version"; \
		awk -v version="$$version" 'BEGIN{p=0} /^\[workspace.package\]$$/{p=1; print; next} /^\[/{p=0} p && /^version = "/{$$0="version = \"" version "\""} {print}' Cargo.toml > Cargo.toml.release.tmp; \
		mv Cargo.toml.release.tmp Cargo.toml; \
		version_changed=1; \
		trap 'status=$$?; if [ $$status -ne 0 ] && [ "$${version_changed:-0}" = 1 ]; then git restore -- Cargo.toml Cargo.lock 2>/dev/null || true; fi; exit $$status' EXIT; \
		echo "Running local musl tests for $$version..."; \
		cargo build -r --target x86_64-unknown-linux-musl; \
		git add Cargo.toml Cargo.lock; \
		git commit -m "Set version to $$version"; \
		version_changed=0; \
		if [ -n "$$(git status --porcelain)" ]; then \
			echo "ERROR: release tests left unexpected repository changes"; \
			git status --short; \
			exit 1; \
		fi; \
		git push origin HEAD; \
		git tag "$$tag"; \
		git push origin "$$tag"; \
	fi; \
	echo "Waiting for GitHub Actions build for $$tag..."; \
	while :; do \
		run_id=$$(gh run list \
			--repo stela2502/lumrik \
			--branch "$$tag" \
			--limit 1 \
			--json databaseId \
			--jq '.[0].databaseId // empty'); \
		[ -n "$$run_id" ] && break; \
		sleep 3; \
	done; \
	conclusion=$$(gh run view "$$run_id" --repo stela2502/lumrik --json conclusion --jq '.conclusion // empty'); \
	if [ "$$conclusion" = "failure" ] || [ "$$conclusion" = "cancelled" ]; then \
		echo "Re-running failed GitHub Actions build $$run_id for $$tag..."; \
		gh run rerun "$$run_id" --repo stela2502/lumrik; \
	fi; \
	gh run watch "$$run_id" \
		--repo stela2502/lumrik \
		--exit-status; \
	image="$(LUMRIK_IMAGE):$$tag"; \
	sif="$(SINGULARITY_LOCAL)/ghcr.io-stela2502-lumrik-$$tag.img"; \
	echo "Pulling $$image"; \
	mkdir -p "$(SINGULARITY_LOCAL)"; \
	rm -f "$$sif"; \
	apptainer pull "$$sif" "docker://$$image"; \
	echo "Copying to COSMOS-SENS cache"; \
	mkdir -p "$(SINGULARITY_SENS)"; \
	cp -v "$$sif" "$(SINGULARITY_SENS)/"; \
	echo "Done: $$tag"

NORN_DIR := $(ROOT)/norn
NORN_CONFIG := $(NORN_DIR)/nextflow.config
NORN_SENS := $(HOME)/sens05_home/norn
.PHONY: lumrik-release-image

.PHONY: norn-deploy-sens
.PHONY: norn_tar_file

norn_tar_file:
	tar -czf norn-source.tar.gz \
		--exclude='target' \
		--exclude='.git' \
		--exclude='*.bam' \
		--exclude='*.bai' \
		--exclude='*.cram' \
		--exclude='*.fastq' \
		--exclude='*.fastq.gz' \
		--exclude='*.fq' \
		--exclude='*.fq.gz' \
		--exclude='*.sif' \
		--exclude='*.rds' \
		--exclude='*.h5' \
		--exclude='*.h5ad' \
		--exclude='*.mtx' \
		--exclude='*.mtx.gz' \
		--exclude='*/whitelists/*.txt.gz' \
		--exclude='*/star_index' \
		--exclude='*/.idea' \
		--exclude='*/.vscode' \
		--exclude='*.swp' \
		norn
	@du -h norn-source.tar.gz



norn-deploy-sens:
	@set -e; \
	echo "=== Norn deployment to COSMOS-SENS ==="; \
	echo ""; \
	\
	cd "$(LUMRIK_DIR)"; \
	tag=$$(git tag --list 'v[0-9]*.[0-9]*.[0-9]*' --sort=-v:refname | head -n1); \
	if [ -z "$$tag" ]; then \
		echo "ERROR: no Lumrik release tag found"; \
		exit 1; \
	fi; \
	echo "Latest Lumrik release: $$tag"; \
	\
	sif="ghcr.io-stela2502-lumrik-$$tag.img"; \
	local_sif="$(SINGULARITY_LOCAL)/$$sif"; \
	\
	if [ ! -f "$$local_sif" ]; then \
		echo "ERROR: matching Lumrik container does not exist:"; \
		echo "  $$local_sif"; \
		echo ""; \
		echo "Run 'make lumrik-release-image' first."; \
		exit 1; \
	fi; \
	\
	cd "$(NORN_DIR)"; \
	\
	current=$$(sed -n \
		"s|.*ghcr.io/stela2502/lumrik:\(v[0-9][0-9.]*\).*|\1|p" \
		"$(NORN_CONFIG)" | head -n1); \
	if [ -z "$$current" ]; then \
		echo "ERROR: could not find Lumrik container version in:"; \
		echo "  $(NORN_CONFIG)"; \
		exit 1; \
	fi; \
	\
	echo "Norn currently uses:    $$current"; \
	\
	bumped=0; \
	if [ "$$current" != "$$tag" ]; then \
		echo "Updating Norn runtime:  $$current -> $$tag"; \
		sed -i \
			"s|ghcr.io/stela2502/lumrik:$$current|ghcr.io/stela2502/lumrik:$$tag|g" \
			"$(NORN_CONFIG)"; \
		bumped=1; \
	else \
		echo "Norn runtime is current."; \
	fi; \
	\
	if ! grep -Fq "ghcr.io/stela2502/lumrik:$$tag" "$(NORN_CONFIG)"; then \
		echo "ERROR: failed to set Norn runtime to Lumrik $$tag"; \
		exit 1; \
	fi; \
	\
	commit_message=""; \
	if [ -n "$$(git status --porcelain)" ]; then \
		echo ""; \
		echo "=== Norn changes ==="; \
		git status --short; \
		echo ""; \
		\
		changed_files=$$(git status --porcelain | awk '{print $$2}'); \
		default_message="Bump Lumrik runtime to $$tag"; \
		\
		if [ "$$changed_files" = "nextflow.config" ] && \
		   grep -Fq "ghcr.io/stela2502/lumrik:$$tag" "$(NORN_CONFIG)"; then \
			commit_message="$$default_message"; \
			echo "Commit message: $$commit_message"; \
		else \
			printf "Commit message [$$default_message]: "; \
			read commit_message; \
			if [ -z "$$commit_message" ]; then \
				commit_message="$$default_message"; \
			fi; \
		fi; \
	else \
		echo ""; \
		echo "Norn repository has no changes to commit."; \
	fi; \
	\
	echo ""; \
	echo "=== Ready ==="; \
	echo "Lumrik release: $$tag"; \
	echo "Lumrik container:"; \
	echo "  $$local_sif"; \
	echo "Norn config:"; \
	echo "  $(NORN_CONFIG)"; \
	if [ -n "$$commit_message" ]; then \
		echo "Norn commit:"; \
		echo "  $$commit_message"; \
	fi; \
	echo ""; \
	printf "Continue? [Y/n] "; \
	read answer; \
	case "$$answer" in \
		n|N|no|NO) echo "Aborted."; exit 1 ;; \
		*) ;; \
	esac; \
	\
	if [ -n "$$commit_message" ]; then \
		echo ""; \
		echo "Committing Norn..."; \
		git add -A; \
		git commit -m "$$commit_message"; \
		echo ""; \
		echo "Pushing Norn..."; \
		git push origin HEAD; \
	fi; \
	\
	echo ""; \
	echo "Deploying to COSMOS-SENS..."; \
	mkdir -p "$(SINGULARITY_SENS)" "$(NORN_SENS)"; \
	rsync -r --progress "$$local_sif" "$(SINGULARITY_SENS)/"; \
	echo "Deploying complete Norn source tree:"; \
	echo "  $(NORN_DIR) -> $(NORN_SENS)"; \
	rsync -r --delete --progress \
		--exclude='.git/' \
		--exclude='.nextflow/' \
		--exclude='.nextflow.log*' \
		--exclude='work/' \
		"$(NORN_DIR)/" "$(NORN_SENS)/"; \
	if [ ! -f "$(NORN_SENS)/main.nf" ] || [ ! -f "$(NORN_SENS)/nextflow.config" ]; then \
		echo "ERROR: incomplete Norn deployment in $(NORN_SENS)"; \
		exit 1; \
	fi; \
	\
	echo ""; \
	echo "Done: Norn + Lumrik $$tag deployed to COSMOS-SENS"


PROJECT ?= lumrik
TARFILE := $(PROJECT)-source.tar.gz

source_tar:
	@test -d "$(PROJECT)" || { echo "ERROR: $(PROJECT) does not exist"; exit 1; }
	tar -czf "$(TARFILE)" \
		--exclude='$(PROJECT)/target' \
		--exclude='$(PROJECT)/.git' \
		--exclude='*.bam' \
		--exclude='*.bai' \
		--exclude='*.cram' \
		--exclude='*.fastq' \
		--exclude='*.fastq.gz' \
		--exclude='*.fq' \
		--exclude='*.fq.gz' \
		--exclude='*.sif' \
		--exclude='*.rds' \
		--exclude='*.h5' \
		--exclude='*.h5ad' \
		--exclude='*.mtx' \
		--exclude='*.mtx.gz' \
		--exclude='*/whitelists/*.txt.gz' \
		--exclude='*/star_index' \
		--exclude='*/.idea' \
		--exclude='*/.vscode' \
		--exclude='*.swp' \
		"$(PROJECT)"
	@du -h "$(TARFILE)"
