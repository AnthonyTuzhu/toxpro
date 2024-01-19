#!/bin/bash
source venv/bin/activate
exec gunicorn --bind 0.0.0.0:443 --worker-tmp-dir /dev/shm --ssl-version TLSv1_2 --certfile /root/ssl/cert.pem --keyfile /etc/letsencrypt/live/toxiverse.com/.pem --workers=2 --timeout 90  --access-logfile - --error-logfile - "app:create_app()"