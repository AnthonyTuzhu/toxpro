#!/bin/bash
source venv/bin/activate
exec gunicorn -b :5000 --access-logfile logs/access.log --error-logfile logs/error.log "app:create_app()"