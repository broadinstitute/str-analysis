TAG1 = weisburd/str-analysis:v0.2.1
TAG2 = weisburd/str-analysis:latest

all: build push

build:
	cd ../ && docker build -t $(TAG1) -t $(TAG2) -f docker/Dockerfile .

push:
	docker push $(TAG1)
	docker push $(TAG2)
