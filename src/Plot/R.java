package Plot;

import AoUtils.AoFile;

import java.util.HashMap;

/**
 * @author AoyueBi
 *
 */
public class R {

    public R(){
        this.changeXaxisLabel();


    }

    /**
     * 在进行箱线图画图的时候，将index具体指的亚种展示出来
     */
    public void changeXaxisLabel(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/WheatVMapII_index_subspecies.txt";
        new AoFile().readheader(infileS);
        HashMap<Integer,String> hm = new AoFile().getHashMapintKey(infileS,0,3);
//        int[] index = {2,3,4,7,10,11,19,13,20,17,14};
        int[] index = {0,10,11,19,13,20,17,14};
        //scale_x_discrete(limits=c("2", "3", "4"),labels = c("wild","Domes","Free"))
        System.out.print("scale_x_discrete(limits=c(");
        for (int i = 0; i < index.length; i++) {
            System.out.print("\"" + index[i] + "\",");
        }
        System.out.print("),labels = c (");
        for (int i = 0; i < index.length; i++) {
            System.out.print("\"" + hm.get(index[i]) + "\",");
        }
        System.out.print("));");


    }

    public static void main(String args[]){
        new R();
    }


}
