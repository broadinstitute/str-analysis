TAG = weisburd/str-analysis-with-eh-and-trgt:latest

all: build push

build:
	docker build -t $(TAG) .

push:
	docker push $(TAG) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
