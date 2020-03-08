from flask import Flask,  jsonify, request, render_template, Blueprint
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
from flask_misaka import Misaka

logging.basicConfig(level=logging.INFO)
logger=logging.getLogger(__name__)
logger.info('Starting wormcells-de...')

flask_app = Flask(__name__)
Misaka(flask_app, math_explicit = True)


tables = Blueprint('tables', __name__, url_prefix='/tables')

# df with the number of cells of each label in each dataset
df = pd.read_csv(flask_app.open_resource('df.csv'))

# to render the table titles better we replace underscores with spaces,
# use non breaking hyphens (&#8209;) and say batch1 instead of just 1
df_nice_names = df.copy()
df_nice_names.columns = df_nice_names.columns.str.replace('_',' ')
df_nice_names.columns = df_nice_names.columns.str.replace('cho-1 1','cho-1 batch1')
df_nice_names.columns = df_nice_names.columns.str.replace('cho-1 2','cho-1 batch2')
df_nice_names.columns = df_nice_names.columns.str.replace('unc-47 2','unc-47 batch2')
df_nice_names.columns = df_nice_names.columns.str.replace('unc-47 1','unc-47 batch1')
df_nice_names.columns = df_nice_names.columns.str.replace('-','&#8209;')

# same for cell type names
# df_nice_names['Cell Type']= df_nice_names['Cell Type'].str.replace('_',' ')


# convert df to dict for sending as json to datatables
dict_df = df_nice_names.to_dict(orient='records')
# convert column names into dict for sending as json to datatables
columns = [{"data": item, "title": item} for item in df_nice_names.columns]

#### datatables ####

@tables.route("/", methods=['GET'])
def clientside_table_content():
    return jsonify({'data': dict_df, 'columns': columns})

flask_app.register_blueprint(tables)

@flask_app.route("/")
def clientside_table():
    return render_template("clientside_table.html")


####

@flask_app.route("/test")
def test():
    return render_template("test.html")

# @flask_app.route("/")
# def index():
#     logger.info('Got a request for index!')
#     return render_template("index.html")

@flask_app.route('/submit', methods=['POST', 'GET'])
def receive_submission():
    logger.info('Got a submission!')
    # answer is a dict of json strings containing selected row and column index numbers
    answer = request.form.to_dict(flat=False)
    print(answer)
    print(df.head())

    #first try is in case submission is from table form
    try:
        # need to convert the json strings to dict, then to a data frame
        # data1 is the selection for the first group, data2 for the second
        data1 = json.loads(answer['data1'][0])
        data1_df = pd.DataFrame.from_dict(data1[0])

        print(data1_df)
        data2 = json.loads(answer['data2'][0])
        data2_df = pd.DataFrame.from_dict(data2[0])

        # now map the index number to experiment name and cell type name
        group1_df = pd.DataFrame()
        group1_df['cell_type1'] = data1_df['row'].map(df['Cell Type'])
        group1_df['experiment1'] = data1_df['column'].map(pd.Series(df.columns.values))
        print(group1_df)

        group2_df = pd.DataFrame()
        group2_df['cell_type2'] = data2_df['row'].map(df['Cell Type'])
        group2_df['experiment2'] = data2_df['column'].map(pd.Series(df.columns.values))
        print(group2_df)

        email = answer['email'][0].strip()
        print(email)

        timestr = time.strftime("%Y%m%d-%H%M%S")
        print(timestr)

        genes = StringIO(json.loads(answer['genes'][0]))
        genes_df = pd.read_csv(genes, names=['selected_genes'])
        print(genes_df)

    # if that doesn't work it's a text submission
    except:
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
    logger.info('the object has been put')
    logger.info(s3filename)
    print()

    ec2 = session.resource('ec2')
    user_data = '''#!/bin/bash
pwd > /home/ubuntu/started_ok.txt
runuser -l  ubuntu -c 'pwd > /home/ubuntu/iamubuntu.txt'   ;
runuser -l  ubuntu -c 'wget -O /home/ubuntu/scvi_de.py https://raw.githubusercontent.com/Munfred/wormcells-de/master/scvi_de.py'   ;
echo 'python3 /home/ubuntu/scvi_de.py ''' + url + ' ' + AWS_S3_ACCESS_KEY + ' ' + AWS_S3_SECRET + ' ' + sendgrid_key + ' ' + sendgrid_name + '''  /home/ubuntu/command.txt'
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

    logger.info('the instance has been created')


    return 'derpderp'

if __name__ == "__main__":
    flask_app.run(host='0.0.0.0')