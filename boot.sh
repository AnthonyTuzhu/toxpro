#!/bin/bash
source venv/bin/activate
exec gunicorn -b :5000 --worker-tmp-dir /dev/shm --workers=2  --access-logfile - --error-logfile - "app:create_app()"