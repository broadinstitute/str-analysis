DOCKERHUB_TAG = weisburd/str-analysis:latest
ARTIFACT_REGISTRY_TAG = us-central1-docker.pkg.dev/cmg-analysis/docker-repo/str-analysis:latest
DOCKERHUB_BASE = $(subst :latest,,$(DOCKERHUB_TAG))
ARTIFACT_REGISTRY_BASE = $(subst :latest,,$(ARTIFACT_REGISTRY_TAG))

all: build push

build:
	docker build -t $(DOCKERHUB_TAG) -t $(ARTIFACT_REGISTRY_TAG) .

push:
	docker push $(DOCKERHUB_TAG) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 | sed 's|^|$(DOCKERHUB_BASE)@|' > sha256_dockerhub.txt
	docker push $(ARTIFACT_REGISTRY_TAG) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 | sed 's|^|$(ARTIFACT_REGISTRY_BASE)@|' > sha256_artifact_registry.txt
