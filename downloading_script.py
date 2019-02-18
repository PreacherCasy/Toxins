from lxml import html
import requests
from bs4 import BeautifulSoup
import pandas as pd
from lxml.etree import fromstring
from Bio import Entrez
from urllib.error import HTTPError
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# taken from , with minor situational additions
class HTMLTableParser:

    def parse_url(self, url):
        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'lxml')
        return [(self.parse_html_table(table)) \
                for table in soup.find_all('table')]

    def parse_html_table(self, table):
        n_columns = 0
        n_rows = 0
        column_names = []
        hyperlinks = []

        for row in table.find_all('tr'):
            td_tags = row.find_all('td')
            if len(td_tags) > 0:
                n_rows += 1
                if n_columns == 0:
                    n_columns = len(td_tags)

            th_tags = row.find_all('th')
            if len(th_tags) > 0 and len(column_names) == 0:
                for th in th_tags:
                    column_names.append(th.get_text())

            hypers = row.find_all('a')
            hyperlinks.append(hypers[0].attrs['href']) if len(hypers) > 0 else hyperlinks.append('NA')

        if len(column_names) > 0 and len(column_names) != n_columns:
            raise Exception("Column titles do not match the number of columns")

        columns = column_names if len(column_names) > 0 else range(0, n_columns)
        df = pd.DataFrame(columns=columns,
                          index=range(0, n_rows), )
        row_marker = 0
        for row in table.find_all('tr'):
            column_marker = 0
            columns = row.find_all('td')
            for column in columns:
                df.iat[row_marker, column_marker] = column.get_text()
                column_marker += 1
            if len(columns) > 0:
                row_marker += 1

        df.columns = list(map(lambda x: x.replace('\n', ''), list(df.iloc[0])))
        df = df.drop(0, axis=0)

        for col in df:
            try:
                df[col] = df[col].astype(float)
            except ValueError:
                pass

        hyperlinks = pd.DataFrame({'NCBI Hyperlink': hyperlinks})
        df = pd.concat([df, hyperlinks.drop(0, axis=0)], axis=1, join='inner')

        return df


init_url = 'http://www.lifesci.sussex.ac.uk/home/Neil_Crickmore/Bt/'
page = requests.get(init_url)
webpage = html.fromstring(page.text)
src_url = webpage.cssselect("frame")[0].attrib["src"]
holo_list = html.fromstring(requests.get(f"{init_url}/{src_url}").text)
tox_list_url = f"{init_url}/{holo_list.xpath('//a/@href')[holo_list.xpath('//a/@href').index('toxins2.html')]}"


hp = HTMLTableParser()
table = hp.parse_url(tox_list_url)[0]
breakdown_list = []

Entrez.email ="271296251017a@gmail.com"

for _, line in table.iterrows():
    ID = line['Name']
    Comment = line['Comment'] if len(line['Comment']) > 1 else None
    Strain = line['Strain / Other ID']
    Accession = line['NCBI Protein'] if len(line['NCBI Protein']) > 1 else line['Acc No.']
    try:
        fetch = Entrez.efetch(db='protein', id=line['Acc No.'], rettype='fasta', retmode='txt')
        read = SeqIO.read(fetch, 'fasta')
        sequence = read.seq if len(set(read.seq.strip())) > 4 else read.seq.translate()
        file = open(f'/home/reverend_casy/Bt_toxins/Downloaded/{ID}.fasta', 'w+')
        SeqIO.write(SeqRecord(sequence, id=f"{ID}, {Comment}", description=f'Strain: {Strain}'), file, format='fasta') if Comment\
        is not None else SeqIO.write(SeqRecord(sequence, id=f"{ID}", description=f"Strain: {Strain}"), file, format='fasta')
        file.close()
        print(f"{ID}: Done")
    except HTTPError:
        try:
            fetch = Entrez.efetch(db='nucleotide', id=Accession, rettype='fasta', retmode='txt')
            read = SeqIO.read(fetch, 'fasta')
            sequence = read.seq if len(set(read.seq.strip())) > 4 else read.seq.translate()
            file = open(f'/home/reverend_casy/Bt_toxins/Downloaded/{ID}.fasta', 'w+')
            SeqIO.write(SeqRecord(sequence, id=f"{ID}, {Comment}", description=f'Strain: {Strain}'), file,
                       format='fasta') if Comment \
                                          is not None else SeqIO.write(
               SeqRecord(sequence, id=f"{ID}", description=f"Strain: {Strain}"), file, format='fasta')
            file.close()
            print(f"{ID}: Done")
        except HTTPError:
            print(f"{ID}: Download Breakdown")
            breakdown_list.append(ID)
            continue

print(breakdown_list)
table.to_csv('~/Bt_toxins/Bt_toxin_table.csv', sep='\t', header=True, index=False)
