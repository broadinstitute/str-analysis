DOCKERHUB_TAG = weisburd/str-analysis:latest
ARTIFACT_REGISTRY_TAG = us-central1-docker.pkg.dev/cmg-analysis/docker-repo/str-analysis:latest

all: build push

build:
	docker build -t $(DOCKERHUB_TAG) -t $(ARTIFACT_REGISTRY_TAG) .

push:
	docker push $(DOCKERHUB_TAG) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256_dockerhub.txt
	docker push $(ARTIFACT_REGISTRY_TAG) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256_artifact_registry.txt
