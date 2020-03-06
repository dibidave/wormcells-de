# wormcells-de



Production deployment is done following this guide from Digital Ocean: 
https://www.digitalocean.com/community/tutorials/how-to-serve-flask-applications-with-gunicorn-and-nginx-on-ubuntu-18-04

It is run with `gunicorn --bind 0.0.0.0:5000 wsgi:flask_app`

`flask_app` is the app name inside `app.py`