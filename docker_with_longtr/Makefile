DOCKERHUB_TAG = weisburd/str-analysis-with-longtr:latest
ARTIFACT_REGISTRY_TAG = us-central1-docker.pkg.dev/cmg-analysis/docker-repo/str-analysis-with-longtr:latest
DOCKERHUB_BASE = $(subst :latest,,$(DOCKERHUB_TAG))
ARTIFACT_REGISTRY_BASE = $(subst :latest,,$(ARTIFACT_REGISTRY_TAG))

all: build push

STR_ANALYSIS_COMMIT ?= main

build:
	docker build --build-arg STR_ANALYSIS_COMMIT=$(STR_ANALYSIS_COMMIT) -t $(DOCKERHUB_TAG) -t $(ARTIFACT_REGISTRY_TAG) .

push:
	docker push $(DOCKERHUB_TAG) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 | sed 's|^|$(DOCKERHUB_BASE)@|' > sha256_dockerhub.txt
	docker push $(ARTIFACT_REGISTRY_TAG) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 | sed 's|^|$(ARTIFACT_REGISTRY_BASE)@|' > sha256_artifact_registry.txt
