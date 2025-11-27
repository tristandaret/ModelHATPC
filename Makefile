SHELL := /bin/bash
.PHONY: build build-dev run compose-up shell clean

IMAGE_NAME := modelhatpc:latest
DEV_IMAGE := modelhatpc:dev

build:
	docker build -t $(IMAGE_NAME) .

build-dev:
	docker build --target dev -t $(DEV_IMAGE) .

build-minimal:
	docker build --target minimal -t $(IMAGE_NAME)-minimal .

run:
	docker run --rm -it -v $(PWD):/app $(IMAGE_NAME) /bin/bash

compose-up:
	docker compose up --build

shell:
	docker run --rm -it -v $(PWD):/app $(IMAGE_NAME) /bin/bash

clean:
	docker image rm -f $(IMAGE_NAME) || true
	docker image rm -f $(DEV_IMAGE) || true
