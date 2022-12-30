#!/bin/bash
source venv/bin/activate
exec rq worker --name toxpro-tasks --url redis://redis:6379/0