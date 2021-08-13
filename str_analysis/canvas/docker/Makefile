TAG1 = weisburd/rfc1-canvas:v0.9
TAG2 = weisburd/rfc1-canvas:latest

all: build push

build:
	cd ../.. && docker build -t $(TAG1) -t $(TAG2) -f canvas/docker/Dockerfile .

push:
	docker push $(TAG1)
	docker push $(TAG2)
