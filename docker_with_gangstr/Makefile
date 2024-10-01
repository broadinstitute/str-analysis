TAG1 = weisburd/str-analysis-with-gangstr:v1.2.10
TAG2 = weisburd/str-analysis-with-gangstr:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
