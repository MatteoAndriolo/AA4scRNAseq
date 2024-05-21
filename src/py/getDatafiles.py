import pathlib
import json

def getfile(name:str="Melanoma",what:str="normed_counts") -> dict[str,str]:
    '''
    This function returns the path to the data files for the given dataset.
    
    Parameters:
    name (str): The name of the dataset. Default is "Melanoma". (options: Melanoma, Exp1, Exp2, Exp3, Myocardial, Mouse)
    what (str): The type of data file to be returned. Default is "normed_counts". (options: all, normed_counts)
    Output:
    dict: A dictionary containing the path to the data files.
    '''
    js=json.load(open(pathlib.Path(__file__).parent.parent/"data"/ "folder_structure.json"))
    if name=="Melanoma":
        print(js["Melanoma"]["normed_counts"])
        return {"normed_counts":pathlib.Path(__file__).parent.parent/"data"/"Melanoma"/js["Melanoma"]["normed_counts"]}
    if name=="Exp1":
        if what=="all":
            out={}
            for k,v in js["AllonKleinLab"]["Experiment1"].items():
                out[k]=pathlib.Path(__file__).parent.parent/"data"/"AllonKleinLab"/"Experiment1"/v
            return out
        else :
            return {"normed_counts": pathlib.Path(__file__).parent.parent/"data"/"AllonKleinLab"/"Experiment1"/js["AllonKleinLab"]["Experiment1"]["normed_counts"]}
    if name=="Exp2":
        if what=="all":
            out={}
            for k,v in js["AllonKleinLab"]["Experiment2"].items():
                out[k]=pathlib.Path(__file__).parent.parent/"data"/"AllonKleinLab"/"Experiment2"/v
            return out
        else :
            return {"normed_counts": pathlib.Path(__file__).parent.parent/"data"/"AllonKleinLab"/"Experiment2"/js["AllonKleinLab"]["Experiment2"]["normed_counts"]}
    if name=="Exp3":
        if what=="all":
            out={}
            for k,v in js["AllonKleinLab"]["Experiment3"].items():
                out[k]=pathlib.Path(__file__).parent.parent/"data"/"AllonKleinLab"/"Experiment3"/v
            return out
        else :
            return {"normed_counts": pathlib.Path(__file__).parent.parent/"data"/"AllonKleinLab"/"Experiment3"/js["AllonKleinLab"]["Experiment3"]["normed_counts"]}
    if name=="Myocardial"|name=="Mouse":
        exit("Not implemented yet")