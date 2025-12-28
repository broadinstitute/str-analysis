#TAG = us-central1-docker.pkg.dev/cmg-analysis/docker-repo/str-analysis:latest
TAG = weisburd/str-analysis:latest


all: build push

build:
	docker build -t $(TAG) .

push:
	docker push $(TAG)
	docker push $(TAG) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
