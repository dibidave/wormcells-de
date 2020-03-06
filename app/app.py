from flask import Flask,  jsonify, request, render_template
import logging
import pandas as pd
import sys
import json
import time
import boto3
import decouple
from io import StringIO
import urllib
from flaskext.markdown import Markdown

logging.basicConfig(level=logging.INFO)
logger=logging.getLogger(__name__)
logger.info('This is a string')


flask_app = Flask(__name__)
Markdown(flask_app)

@flask_app.route("/test")
def test():
    return render_template("test.html")

@flask_app.route("/")
def index():
    logger.info('Got a request for index!')
    return render_template("index.html")

@flask_app.route('/submit', methods=['POST', 'GET'])
def receive_submission():
    logger.info('Got a submission!')
    # answer is a dict of json strings containing selected row and column index numbers
    answer = request.form.to_dict(flat=False)
    print(answer)

    # need to convert the json strings to stringio, then to a data frame
    # data1 is the selection for the first group, data2 for the second

    # print(email)

    data1 = StringIO(json.loads(answer['data1'][0]))
    # print(data1)
    group1_df = pd.read_csv(data1, names=['cell_type1', 'experiment1'])
    # print(data1_df)

    data2 = StringIO(json.loads(answer['data2'][0]))
    group2_df = pd.read_csv(data2, names=['cell_type2', 'experiment2'])

    genes = StringIO(json.loads(answer['genes'][0]))
    genes_df = pd.read_csv(genes, names=['selected_genes'])
    print(genes_df)

    # now map the index number to experiment name and cell type name


    email = answer['email'][0].strip()
    print(email)

    timestr = time.strftime("%Y%m%d-%H%M%S")
    print(timestr)

    s3filename = 'submissions/' + email + '%' + timestr + '.csv'

    selected_groups_df = pd.concat([group1_df, group2_df, genes_df], axis=1)
    print(selected_groups_df)
    AWS_S3_ACCESS_KEY = decouple.config('AWS_S3_ACCESS_KEY')
    AWS_S3_SECRET = decouple.config('AWS_S3_SECRET')
    sendgrid_key = decouple.config('sendgrid_key')
    sendgrid_name = decouple.config('sendgrid_name')

    # print(AWS_S3_SECRET)
    # prinzt(AWS_S3_ACCESS_KEY)
    csv_buffer = StringIO()
    selected_groups_df.to_csv(csv_buffer)

    client = boto3.client('s3',
                          aws_access_key_id=AWS_S3_ACCESS_KEY,
                          aws_secret_access_key=AWS_S3_SECRET,
                          region_name='us-east-2'
                          )
    client.put_object(
        Body=csv_buffer.getvalue(),
        Bucket='scvi-differential-expression',
        Key=s3filename,
        ACL='public-read'
    )
    session = boto3.Session(region_name='us-east-2',
                            aws_access_key_id=AWS_S3_ACCESS_KEY,
                            aws_secret_access_key=AWS_S3_SECRET)

    # s3 = session.resource('s3')
    # bucket = s3.Bucket('scvi-differential-expression')
    # bucket.Acl().put(ACL='public-read')
    # bucket.upload_file(csv_buffer.getvalue(), s3filename)


    url = 'https://scvi-differential-expression.s3.us-east-2.amazonaws.com/' + urllib.parse.quote(s3filename)
    print('the objeoct has been put')
    print(s3filename)
    print()

    ec2 = session.resource('ec2')
    user_data = '''#!/bin/bash
pwd > /home/ubuntu/started_ok.txt
runuser -l  ubuntu -c 'pwd > /home/ubuntu/iamubuntu.txt'   ;
runuser -l  ubuntu -c 'wget -O /home/ubuntu/scvi_de.py https://raw.githubusercontent.com/Munfred/wormcells-de/master/scvi_de.py'   ;
echo 'python3 /home/ubuntu/scvi_de.py ''' + url + ' ' + AWS_S3_ACCESS_KEY + ' ' + AWS_S3_SECRET + ' ' + sendgrid_key + ' ' + sendgrid_name + '''  home/ubuntu/command.txt'
runuser -l  ubuntu -c 'python3 /home/ubuntu/scvi_de.py ''' + url + ' ' + AWS_S3_ACCESS_KEY + ' ' + AWS_S3_SECRET + ' ' + sendgrid_key + ' ' + sendgrid_name + ''' ;'
echo "sudo halt" 

'''

    print(user_data)

    # create a new EC2 instance
    instances = ec2.create_instances(
        ImageId='ami-032240eb155129553',
        MinCount=1,
        MaxCount=1,
        InstanceType='r5d.4xlarge',
        UserData=user_data,
        KeyName='ec2-keypair'
    )

    print('the instance has been created')


    return 'derpderp'

if __name__ == "__main__":
    app.run(host='0.0.0.0')