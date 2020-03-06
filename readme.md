# wormcells-de



Production deployment is done following this guide from Digital Ocean: 
https://www.digitalocean.com/community/tutorials/how-to-serve-flask-applications-with-gunicorn-and-nginx-on-ubuntu-18-04

It is run with `gunicorn --bind 0.0.0.0:5000 wsgi:flask_app`

`flask_app` is the app name inside `app.py`

## Setting up a new server instance with ubuntu1 18.04

```
sudo apt update
sudo apt install python3-pip python3-dev build-essential libssl-dev libffi-dev python3-setuptools
sudo apt install python3-venv

cd ~
git clone https://github.com/Munfred/wormcells-de.git
cd wormcells-de
python3.6 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install wheel
pip install gunicorn flask

gunicorn --bind 0.0.0.0:5000 wsgi:flask_app
deactivate
```

### Configuring service
Then, following the the tutorial, create service file 
```
sudo nano /etc/systemd/system/wormcells.service
```
Contents:

```
[Unit]
Description=Gunicorn instance to serve wormcells
After=network.target

[Service]
User=ubuntu
Group=www-data
WorkingDirectory=/home/ubuntu/wormcells-de/
Environment="PATH=/home/ubuntu/wormcells-de/venv/bin"
ExecStart=/home/ubuntu/wormcells-de/venv/bin/gunicorn --workers 5 --bind unix:wormcells.sock -m 007 wsgi:flask_app

[Install]
WantedBy=multi-user.target
```

Then
```
sudo systemctl start wormcells
sudo systemctl enable wormcells
```

### Configuring Nginx

Install Nginx, create server block configuration file:
```
sudo apt install nginx

sudo nano /etc/nginx/sites-available/wormcells
```

Contents:
```
server {
    listen 80;
    server_name wormcells de.wormcells.com;

    location / {
        include proxy_params;
        proxy_pass http://unix:/home/ubuntu/wormcells-de/wormcells.sock;
    }
}
```


To enable the Nginx server block configuration link the file to the sites-enabled directory:

```
sudo ln -s /etc/nginx/sites-available/wormcells /etc/nginx/sites-enabled
```
With the file in that directory, you can test for syntax errors:
```
sudo nginx -t
```

If this returns without indicating any issues, restart the Nginx process to read the new configuration:

** MAKE SURE NOTHING ELSE IS RUNNING ON PORT 80 **
```
sudo systemctl restart nginx
```

