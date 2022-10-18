TAG1 = weisburd/str-analysis:v0.9.8
TAG2 = weisburd/str-analysis:latest

all: build push

build:
	cd ../ && docker build -t $(TAG1) -t $(TAG2) -f docker/Dockerfile .

push:
	docker push $(TAG1)
	docker push $(TAG2) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
