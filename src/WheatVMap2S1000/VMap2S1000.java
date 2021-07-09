package WheatVMap2S1000;

import WheatGeneticLoad.VariantsSum;

public class VMap2S1000 {
    public VMap2S1000(){
        this.snpAnnotationBuild();


    }

    public void snpAnnotationBuild(){
        new VariantsSum().extractInfoFromVMap2();
        new VariantsSum().mkExonVCF();
        new VariantsSum().mkExonAnnotation2();
        new VariantsSum().addSift();
        new VariantsSum().addAncestral();
        new VariantsSum().addDerived_SIFT();
        new VariantsSum().addDAF();
        new VariantsSum().addGerp();
        new VariantsSum().mergeExonSNPAnnotation();


    }


}
