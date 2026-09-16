DOCKERHUB_TAG = weisburd/str-analysis:latest
DOCKERHUB_BASE = $(subst :latest,,$(DOCKERHUB_TAG))

all: build push

STR_ANALYSIS_COMMIT ?= main

build:
	docker build --build-arg STR_ANALYSIS_COMMIT=$(STR_ANALYSIS_COMMIT) -t $(DOCKERHUB_TAG) .

push:
	docker push $(DOCKERHUB_TAG) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 | sed 's|^|$(DOCKERHUB_BASE)@|' > sha256_dockerhub.txt
