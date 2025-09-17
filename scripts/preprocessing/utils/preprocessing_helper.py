import yaml
import pandas as pd

def load_config(config_file_path):
    try:
        with open(config_file_path, 'r') as stream:
            config = yaml.safe_load(stream)
        return config
    except yaml.YAMLError as exc:
        print(f"Error loading YAML file: {exc}")
        return None
    
def prep_PROs_for_mplus(config_data):
    # Import COMFORT questionnaire data
    path = config_data['default']['data_raw'] + config_data['default']['file_PRO']
    HADS_sheet = "HADS"
    df_hads = pd.read_excel(io=path, sheet_name=HADS_sheet)
    df_hads = df_hads[['Participant #', 'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13', 'E14']]
    df_hads.rename(columns={'Participant #': 'id'}, inplace=True)

    SAGIS_sheet = "SAGIS"
    df_sagis = pd.read_excel(io=path, sheet_name=SAGIS_sheet)
    df_sagis.drop(['Participant #.1', 'CASE/Control', 'CASE/CONTROL', 'Epigastric', 'IBS-D', 'Acid', 'Vomiting', 'Constipation'], axis=1, inplace=True)
    dict_SAGIS = {('1. Belching with acid taste/heartburn/burning sensation in the oesophagus '
                '(tube joining mouth and stomach)'): 'sagis_1',
                '2. Dysphagia (difficulty swallowing)': 'sagis_2',
                '3. Fullness (feeling of congestion of food without relation to prior food intake)': 'sagis_3',
                ('4. Early satiety (stomach is overfilled soon after starting to eat, disproportional to the '
                'quantity of food taken, so that food cannot be finished)'): 'sagis_4',
                '5. Postprandial pain or discomfort (upper abdominal symptoms start or get worse after meals)': 'sagis_5',
                ('6. Epigastric pain/upper abdominal pain (pain between the belly button and chest/pain noticeable '
                'in the upper abdomen)'): 'sagis_6',
                '7. Retrosternal discomfort (unpleasant feeling behind the middle of the chest, painful or drawing)':
                    'sagis_7',
                '8. Pain or discomfort prior to bowel movement': 'sagis_8',
                '9. Difficulty with emptying the bowel': 'sagis_9',
                '10. Constipation': 'sagis_10',
                '11. Loose stools': 'sagis_11',
                '12. Incontinence': 'sagis_12',
                '13. Urgency to empty the bowel': 'sagis_13',
                '14. Diarrhoea': 'sagis_14',
                '15. Loss of appetite (listless for food intake)': 'sagis_15',
                '16. Abdominal cramps (spasmodic or colic like stomach pain without specified localisation)': 'sagis_16',
                '17. Sickness (discomfort combined with the impression for the need to vomit)': 'sagis_17',
                '18. Nausea (urgent feeling of the need to vomit)': 'sagis_18',
                '19. Vomiting (vomiting of mucus and gastric contents, or strong unproductive retching)': 'sagis_19',
                '20. Bloating (feeling of distension and excessive gas in the abdomen)': 'sagis_20',
                '21. Excessive gas and passing of wind': 'sagis_21',
                '22. Excessive belching': 'sagis_22',
                'Participant #': 'id'}
    df_sagis.rename(columns=dict_SAGIS, inplace=True)

    PROMIS_sheet = "PROMIS"
    df_promis = pd.read_excel(io=path, sheet_name=PROMIS_sheet)
    df_promis.drop(['CASE/CONTROL', 'CASE/Control', 'GENDER', 'AGE'], axis=1, inplace=True)
    dict_PROMIS = {'Participant #': 'id',
                'P_Anxiety': 'pr_anx',
                'P_Bloating': 'pr_blo',
                'P_Constipation': 'pr_con',
                'P_DS': 'pr_ds',
                'P_Depression': 'pr_dep',
                'P_Diarrhoea': 'pr_dia',
                'P_Pain': 'pr_pain',
                'P_Reflux': 'pr_ref'}
    df_promis.rename(columns=dict_PROMIS, inplace=True)

    # Combine all data
    df = df_hads.merge(df_promis, how='outer', on='id')
    df = df.merge(df_sagis, how='outer', on='id')
                
    # Replace missing values with -999
    df = df.fillna(-999)
    return df


def fig_1_preprocessing(config_data):
    df_amino = pd.read_csv(config_data['default']['data_processed']+"comfort_amino_acids_cluster.csv", index_col = 0).T
    df_bile = pd.read_csv(config_data['default']['data_processed']+"comfort_bile_acids_cluster.csv", index_col = 0).T
    df_f_metab = pd.read_csv(config_data['default']['data_processed']+"comfort_fmetab_cluster.csv", index_col = 0).T
    df_gene = pd.read_csv(config_data['default']['data_processed']+"comfort_gene_cluster.csv", index_col = 0).T
    df_taxa = pd.read_csv(config_data['default']['data_processed']+"comfort_taxa_cluster.csv", index_col = 0).T
    df_scfa = pd.read_csv(config_data['default']['data_processed']+"comfort_organic_acids_cluster.csv", index_col = 0).T
    df_p_metab = pd.read_csv(config_data['default']['data_processed']+"comfort_pmetab_cluster.csv", index_col = 0).T

    df_hads = pd.read_excel(io=config_data['default']['data_raw']+config_data['default']['file_PRO'], sheet_name="HADS")
    df_hads['Participant'] = df_hads.loc[:, 'Participant #']
    df_hads_data = df_hads.iloc[:, 5:19]
    df_hads_data.index = df_hads['Participant']
    df_hads_data.dropna(axis=0, how='all', inplace=True)
    df_hads_data = df_hads_data.notnull().astype("int")
    hads_present = df_hads_data.all(axis=1)
    hads_present = hads_present.astype(int)
    hads_present.replace(0,0.7, inplace=True)
    df_hads_present = pd.concat([df_hads_data.index.to_series(), hads_present], axis=1)
    df_hads_present.rename(columns={0: 'HADS', 'Participant': 'Participants'}, inplace=True)

    df_sagis = pd.read_excel(io=config_data['default']['data_raw']+config_data['default']['file_PRO'], sheet_name="SAGIS")
    df_sagis['Participant'] = df_sagis.loc[:, 'Participant #']
    df_sagis_data = df_sagis.iloc[:, 4:25]
    df_sagis_data.index = df_sagis['Participant']
    df_sagis_data.dropna(axis=0, how='all', inplace=True)
    df_sagis_data = df_sagis_data.notnull().astype("int")
    sagis_present = df_sagis_data.all(axis=1)
    sagis_present = sagis_present.astype(int)
    sagis_present.replace(0,0.7, inplace=True)
    df_sagis_present = pd.concat([df_sagis_data.index.to_series(), sagis_present], axis=1)
    df_sagis_present.rename(columns={0: 'SAGIS', 'Participant': 'Participants'}, inplace=True)

    df_promis = pd.read_excel(io=config_data['default']['data_raw']+config_data['default']['file_PRO'], sheet_name="PROMIS")
    df_promis['Participant'] = df_promis.loc[:, 'Participant #']
    df_promis_data = df_promis.iloc[:, 5:13]
    df_promis_data.index = df_promis['Participant']
    df_promis_data.dropna(axis=0, how='all', inplace=True)
    df_promis_data = df_promis_data.notnull().astype("int")
    promis_present = df_promis_data.all(axis=1)
    promis_present = promis_present.astype(int)
    promis_present.replace(0,0.7, inplace=True)
    df_promis_present = pd.concat([df_promis_data.index.to_series(), promis_present], axis=1)
    df_promis_present.rename(columns={0: 'PROMIS', 'Participant': 'Participants'}, inplace=True)
    df_promis_present.reset_index()

    df_scfa = df_scfa.notnull().astype("int")
    scfa_present = df_scfa.all(axis=1)
    scfa_present = scfa_present.astype(int)
    scfa_present.replace(0,0.7, inplace=True)
    df_scfa_present = pd.concat([df_scfa.index.to_series(), scfa_present], axis=1)
    df_scfa_present.rename(columns={1: 'SCFA', 0: 'Participants'}, inplace=True)
    df_scfa_present['Participants'] = df_scfa_present['Participants'].astype(int)

    df_bile = df_bile.notnull().astype("int") 
    bile_present = df_bile.all(axis=1)
    bile_present = bile_present.astype(int)
    bile_present.replace(0,0.7, inplace=True)
    df_bile_present = pd.concat([df_bile.index.to_series(), bile_present], axis=1)
    df_bile_present.rename(columns={1: 'bile', 0: 'Participants'}, inplace=True)
    df_bile_present['Participants'] = df_bile_present['Participants'].astype(int)

    df_amino = df_amino.notnull().astype("int") 
    amino_present = df_amino.all(axis=1)
    amino_present = amino_present.astype(int)
    amino_present.replace(0,0.7, inplace=True)
    df_amino_present = pd.concat([df_amino.index.to_series(), amino_present], axis=1)
    df_amino_present.rename(columns={1: 'amino', 0: 'Participants'}, inplace=True)
    df_amino_present['Participants'] = df_amino_present['Participants'].astype(int)

    df_fmetab = df_f_metab.notnull().astype("int") 
    fmetab_present = df_fmetab.all(axis=1)
    fmetab_present = fmetab_present.astype(int)
    fmetab_present.replace(0,0.7, inplace=True)
    df_fmetab_present = pd.concat([df_fmetab.index.to_series(), fmetab_present], axis=1)
    df_fmetab_present.rename(columns={1: 'fmetab', 0: 'Participants'}, inplace=True)
    df_fmetab_present['Participants'] = df_fmetab_present['Participants'].astype(int)

    df_pmetab = df_p_metab.notnull().astype("int")
    pmetab_present = df_pmetab.all(axis=1)
    pmetab_present = pmetab_present.astype(int)
    pmetab_present.replace(0,0.7, inplace=True)
    df_pmetab_present = pd.concat([df_pmetab.index.to_series(), pmetab_present], axis=1)
    df_pmetab_present.rename(columns={1: 'pmetab', 0: 'Participants'}, inplace=True)
    df_pmetab_present['Participants'] = df_pmetab_present['Participants'].astype(int)

    df_gene = df_gene.notnull().astype("int")
    gene_present = df_gene.all(axis=1)
    gene_present = gene_present.astype(int)
    gene_present.replace(0,0.7, inplace=True)
    df_gene_present = pd.concat([df_gene.index.to_series(), gene_present], axis=1)
    df_gene_present.rename(columns={1: 'gene', 0: 'Participants'}, inplace=True)
    df_gene_present['Participants'] = df_gene_present['Participants'].astype(int)    

    df_taxa = df_taxa.notnull().astype("int")
    taxa_present = df_taxa.all(axis=1)
    taxa_present = taxa_present.astype(int)
    taxa_present.replace(0,0.7, inplace=True)
    df_taxa_present = pd.concat([df_taxa.index.to_series(), taxa_present], axis=1)
    df_taxa_present.rename(columns={1: 'taxa', 0: 'Participants'}, inplace=True)
    df_taxa_present['Participants'] = df_taxa_present['Participants'].astype(int)

    # Combine all data indicators
    df_present = pd.merge(df_bile_present, df_scfa_present, how="outer", on="Participants")
    df_present = pd.merge(df_present, df_amino_present, how="outer", on="Participants")
    df_present = pd.merge(df_present, df_pmetab_present, how="outer", on="Participants")
    df_present = pd.merge(df_present, df_fmetab_present, how="outer", on="Participants")
    df_present = pd.merge(df_present, df_taxa_present, how="outer", on="Participants")
    df_present = pd.merge(df_present, df_gene_present, how="outer", on="Participants")
    df_present = pd.merge(df_present, df_sagis_present, how="outer", on="Participants")
    df_present = pd.merge(df_present, df_promis_present, how="outer", on="Participants")
    df_present = pd.merge(df_present, df_hads_present, how="outer", on="Participants")

    # Sort by participants
    df_present.sort_values(by=['Participants'], inplace=True)

    # Transpose matrix and make participant number headers
    df_present = df_present.transpose()
    df_present.columns = df_present.iloc[0]
    df_present = df_present[1:]

    df_present.fillna(0, inplace=True)

    return df_present