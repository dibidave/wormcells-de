import os
import numpy as np
import pandas as pd
from scvi.dataset import GeneExpressionDataset
from scvi.models import VAE
from scvi.inference import UnsupervisedTrainer
import torch
import anndata
import plotly.graph_objects as go
import io
import requests
import urllib
from io import StringIO
import boto3
import datetime
import sys
from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail

import datetime
dt_started = datetime.datetime.utcnow()

url = (sys.argv[1])
print(datetime.datetime.now())
print('	### ### ###  starting to process...')
print(url)

AWS_S3_ACCESS_KEY = (sys.argv[2])

AWS_S3_SECRET = (sys.argv[3])

sendgrid_key = (sys.argv[4])

sendgrid_name = (sys.argv[5])

save_path = "/home/ubuntu/"

sys.stdout = open('log_file.txt', 'w')
sys.stderr = sys.stdout
# if True:
try:
    vae_file_name = 'quis9h95ckgwf0vo8ozsxif72qxsivfi.pkl'
    adata = anndata.read('./n5xlzu8zum9c0nq3blinxg5ml9s6nab5.h5ad')
#     adata.var = adata.var.reset_index()
#     adata.var.index = adata.var['gene_id']
#     adata.var.index = adata.var.index.rename('index')

    gene_dataset = GeneExpressionDataset()

    # we provide the `batch_indices` so that scvi can perform batch correction
    gene_dataset.populate_from_data(
        adata.X,
        gene_names=adata.var.index.values,
        cell_types=adata.obs['cell_type_infection'].values,
        batch_indices=adata.obs['sample'].cat.codes.values,
    )

    vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches)

    trainer = UnsupervisedTrainer(
        vae,
        gene_dataset,
        train_size=0.75,
        use_cuda=False,
        frequency=1,
    )

    # LOAD
    full_file_save_path = os.path.join(save_path, vae_file_name)
    trainer.model.load_state_dict(torch.load(full_file_save_path))
    trainer.model.eval()
    print('	### ### ###  loaded vae')
    print(datetime.datetime.now())

    # n_epochs = 5
    # lr = 0.001
    # full_file_save_path = os.path.join(save_path, vae_file_name)
    # trainer.train(n_epochs=n_epochs, lr=lr)
    # torch.save(trainer.model.state_dict(), full_file_save_path)
    # train_test_results = pd.DataFrame(trainer.history).rename(columns={'elbo_train_set':'Train', 'elbo_test_set':'Test'})
    # print(train_test_results)

    full = trainer.create_posterior(trainer.model, gene_dataset, indices=np.arange(len(gene_dataset)))
    latent, batch_indices, labels = full.sequential().get_latent()
    batch_indices = batch_indices.ravel()

    print('	### ### ###  computed full posterior')


    print('	### ### ###  url = ', url)
    # read submission csv and fetch selected cells
    submission = pd.read_csv(io.StringIO(requests.get(url).content.decode('utf-8')), index_col=0)
    selected_cells_csv_string = submission.to_csv(index=False).replace('\n', '<br>')

    # reconstruct user email from submission url
    email = url.split('https://aavcells-de.s3.us-west-2.amazonaws.com/submissions/')[1]
    email = email.split('%25')[0]
    email = email.replace('%40', '@')

    # create masks for cell selection according to submission
    # start with all entries false
    cell_idx1 = adata.obs['cell_type_infection'] == '000000'

    for _, entry in submission[['cell_type1', 'experiment1']].dropna().iterrows():
        cell = entry['cell_type1'].strip()
        experiment = entry['experiment1'].strip()
        curr_boolean = (adata.obs['cell_type_infection'] == cell) & (adata.obs['sample'] == experiment)
        cell_idx1 = (cell_idx1 | curr_boolean)

    cell_idx2 = adata.obs['cell_type_infection'] == '000000'

    for _, entry in submission[['cell_type2', 'experiment2']].dropna().iterrows():
        cell = entry['cell_type2'].strip()
        experiment = entry['experiment2'].strip()
        curr_boolean = (adata.obs['cell_type_infection'] == cell) & (adata.obs['sample'] == experiment)
        cell_idx2 = (cell_idx2 | curr_boolean)

    print('	### ### ###  computed idxs')

    n_samples = 10000
    M_permutation = 10000

    de_change = full.differential_expression_score(
        idx1=cell_idx1.values,  # we use the same cells as chosen before
        idx2=cell_idx2.values,
        mode='change',  # set to the new change mode
        n_samples=n_samples,
        M_permutation=M_permutation,
    )

    print('	### ### ###  finished DE!')
    print(datetime.datetime.now())

    # manipulate the DE results for plotting

    # we use the `mean` entru in de_chage, it is the scVI posterior log2 fold change
    de_change['log10_pvalue'] = np.log10(de_change['proba_not_de'])

    # we take absolute values of the first bayes factor as the one to use on the volcano plot
    # bayes1 and bayes2 should be roughtly the same, except with opposite signs
    de_change['abs_bayes_factor'] = np.abs(de_change['bayes_factor'])
    de_change = de_change.join(adata.var, how='inner')

    # manipulate the DE results for plotting
    de = de_change.copy()

    de['gene_description'] = ''
    de['gene_name'] = de.index
    de['gene_id'] = de.index
    
    de['gene_color'] = 'rgba(100, 100, 100, 0.2)'
    for gene in submission['selected_genes'].dropna().values:
        gene = gene.strip()
        de['gene_color'][de['gene_name'].str.contains(gene)] = 'rgba(0, 0,255, 1)'
        de['gene_color'][de['gene_id'].str.contains(gene)] = 'rgba(0, 0,255, 1)'
    
    de_result_csv = de[
        ['proba_not_de', 'log10_pvalue', 'bayes_factor', 'lfc_mean', 'lfc_median', 'lfc_std', 'lfc_min', 'lfc_max',
         'gene_id', 'gene_name', 'gene_description']]

    print('	### ### ###  Creating plot')

    try:
        jobname = submission['job_name'][0]
    except:
        jobname = ' '

    de['gene_description_html'] = de.index
    string_bf_list = [str(bf) for bf in np.round(de['bayes_factor'].values, 3)]
    de['bayes_factor_string'] = string_bf_list

    fig = go.Figure(
        data=go.Scatter(
            x=de["lfc_mean"].round(3)
            , y=-de["log10_pvalue"].round(3)
            , mode='markers'
            , marker=dict(color=de['gene_color'])
            , hoverinfo='text'
            , text=de['gene_description_html']
            , customdata=de.gene_name.astype(str) + '<br>' + de.gene_id.values + \
                         '<br>Bayes Factor: \t' + de.bayes_factor_string + \
                         '<br>-log10(p-value): \t' + de["log10_pvalue"].round(3).astype(str) + \
                         '<br>log2 FC mean: \t' + de["lfc_mean"].round(3).astype(str) + \
                         '<br>log2 FC median: \t' + de["lfc_median"].round(3).astype(str) + \
                         '<br>log2 FC std: \t' + de["lfc_std"].round(3).astype(str)
            , hovertemplate='%{customdata} <br><extra>%{text}</extra>'
        )
        , layout={
            "title": {"text":
                          "Differential expression on AAV single cell data  <br> " + jobname + " <br> <a href=" + url + ">" + url + '</a>'
                , 'x': 0.5
                      }
            , 'xaxis': {'title': {"text": "log2 fold change"}}
            , 'yaxis': {'title': {"text": "-log10(p-value)"}}
        }
    )

    ### SAVE RESULTS AND SEND MAIL ####

    # construct the filename for saving the results
    filename = url.split('https://aavcells-de.s3.us-west-2.amazonaws.com/submissions/')[1]
    filename = filename.replace('%40', '@')
    filename = filename.replace('%25', '@')
    filename = filename.replace('.csv', '')

    csv_buffer = StringIO()
    de_result_csv.to_csv(csv_buffer)
    csvfilename = 'csv/' + filename + '-results.csv'
    htmlfilename = 'plots/' + filename + '-results.html'
    print('	### ### ###  Putting files in s3...')

    client = boto3.client('s3')

    client.put_object(
        Body=csv_buffer.getvalue(),
        Bucket='aavcells-de',
        Key=csvfilename,
        ACL='public-read'
    )

    html_buffer = StringIO()
    fig.write_html(html_buffer, auto_open=True)

    client.put_object(
        Body=html_buffer.getvalue(),
        Bucket='aavcells-de',
        Key=htmlfilename,
        ACL='public-read'
    )

    csv_url = 'https://aavcells-de.s3.us-west-2.amazonaws.com/' + urllib.parse.quote(csvfilename)
    html_url = 'https://aavcells-de.s3.us-west-2.amazonaws.com/' + urllib.parse.quote(htmlfilename)
    print('	### ### ###  Files uploaded successfully')
    dt_ended = datetime.datetime.utcnow()
    total_time = str(int((dt_ended - dt_started).total_seconds()))

    email_body = f' Your aavcells-de single cell differential expression results: <br> {jobname} <br><br> <a href="{csv_url}">CSV file with results</a>  <br> <a href="{html_url}">Vocano plot</a>  <br> <a href="{url}">Your cell selection</a> <br> Processing time: {total_time}s <br> <br> <br> Thanks <br> Eduardo and David'
    print('	### ### ###  Email created')

    message = Mail(
        from_email='dibidave@caltech.edu',
        to_emails=email,
        subject='AAV single-cell differential expression results',
        html_content=email_body)
    sg = SendGridAPIClient(sendgrid_key)
    response = sg.send(message)
    print(response.status_code)
    print(response.body)
    print(response.headers)
    print('	********************* Email sent ********************')
    print('	****************** Total time (s):', end = ' ')
    print(total_time)

    print('Terminating... ')
    instance_id = requests.get("http://169.254.169.254/latest/meta-data/instance-id").text

    session = boto3.Session(region_name='us-west-2')
    ec2 = session.resource('ec2', region_name='us-west-2')
    ec2.instances.filter(InstanceIds=[instance_id]).terminate()
# if False:
except:
    print('	XXXXXXXXXXXXXXXX SOMETHING FAILED XXXXXXXXXXXXXXX')
    e = sys.exc_info()[0]
    print(e)
    filename = url.split('https://aavcells-de.s3.us-west-2.amazonaws.com/submissions/')[1]
    filename = filename.replace('%40', '@')
    filename = filename.replace('%25', '@')
    filename = filename.replace('.csv', '')

    logs = open('log_file.txt', 'r').read()

    logs_buffer = StringIO()
    logs_buffer.write(logs)

    csv_buffer = StringIO()
    logsfilename = 'logs/' + filename + '-logs.txt'
    print('	### ### ###  Putting log in s3...')

    client = boto3.client('s3')

    client.put_object(
        Body=logs_buffer.getvalue(),
        Bucket='aavcells-de',
        Key=logsfilename,
        ACL='public-read'
    )
    print('	### ### ###  Log in s3. Sending email...')

    message = Mail(
        from_email='dibidave@caltech.edu',
        to_emails='dibidave@caltech.edu',
        subject='SOMETHING WENT WRONG!!!!111',
        html_content=logs)

    sg = SendGridAPIClient(sendgrid_key)
    response = sg.send(message)
    print(response.status_code)
    print(response.body)
    print(response.headers)

    print('	### ### ###  Email sent. Done.')

    print('Terminating... ')
    instance_id = requests.get("http://169.254.169.254/latest/meta-data/instance-id").text
    
    session = boto3.Session(region_name='us-west-2')
    ec2 = session.resource('ec2', region_name='us-west-2')
    ec2.instances.filter(InstanceIds = [instance_id]).terminate()
