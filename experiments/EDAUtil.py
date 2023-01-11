import os

class EDAUtil(object):
    def __init__(self):
        pass
    
    def get_rbps(self, directory):
        rbps = []
        for file in os.listdir(f'{directory}'):
            with open(f'{directory}/{file}') as embeddings:
                for embedding in embeddings:
                    rbp_id = embedding.split(',')[0]
                    if rbp_id != 'ID':
                        rbps.append(rbp_id)
                
        return rbps