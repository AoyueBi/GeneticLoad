package AoUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class AoColor {

    public AoColor(){
        this.getCol();

    }


    /**
     * e.g.
     */
    public void getCol(){
//        String[] in = {"AB","ABD","D"};
//        AoColor.genomeType(in);

//        String[] in = {"Ae.tauschii","Cultivar","Domesticated_emmer","Free_threshing_tetraploid","Landrace","OtherHexaploid","OtherTetraploid","Wild_emmer"};
//        String[] in = {"Cultivar","Domesticated_emmer","Free_threshing_tetraploid","Landrace","OtherHexaploid","OtherTetraploid","Wild_emmer"};
//        String[] in = {"Ae.tauschii","Cultivar","Landrace","OtherHexaploid"};
//        String[] in = {"Domesticated_emmer","Free_threshing_tetraploid","OtherTetraploid","Wild_emmer"};

//        AoColor.subspecies(in);

//        String[] in = {"Africa" , "America" ,"Asia" ,   "Europe" , "Oceania"};
//        String[] in = {"Africa",  "America" ,"Asia"   , "Europe"};
        String[] in = {"America" ,"Asia" ,   "Europe" };

        AoColor.continent(in);

    }


    //分组2为：大洋洲 非洲 北美洲 南美洲 欧洲 亚洲
    //"Oceania","Africa","North America","South America","Europe","Asia"
    //颜色为："#F1E1FF","#F4D03F","#5DADE2","#7B241C","#FF9900","#82C782"
    //
    /**
     *
     * @param input
     * @return
     */

    public static String continent (String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#7dbde8", "#7B241C", "#dbb3ff","#FF9900", "#82C782"};
        String[] key = {"Oceania","Africa","America","Europe","Asia"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }

    /**
     *
     * @param input
     * @return
     */

    public static String subspecies (String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#87cef9","#9900ff","#7f5701","#016699","#fc6e6e","#fe63c2","#00f3ff","#ffd702"};
        String[] key = {"Ae.tauschii","Cultivar","Domesticated_emmer","Free_threshing_tetraploid","Landrace","OtherHexaploid","OtherTetraploid","Wild_emmer"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }


    /**
     *
     * @param input
     * @return
     */

    public static String genoType2(String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#ffd702","#fc6e6e","#87cef9"};
        String[] key = {"AABB","AABBDD","DD"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }

    public static String genomeType(String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#ffd702","#fc6e6e","#87cef9"};
        String[] key = {"AB","ABD","D"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }

    public static String models7(String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#fd8582","#ffcbcc","#967bce","#cccdfe","darkgray","#4bcdc6","#cdfffc"};
        String[] key = {"M100","M011","M010","M101","M111","M001","M110"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }

    /**
     * 将因子变量进行自然顺序的排序，然后再调用颜色值
     * @param input
     */
    public static void sort(String[] input){
//        String[] order = {"M001","M010","M011", "M100", "M101", "M110", "M111"};
        List<String> l = new ArrayList<>();
        for (int i = 0; i < input.length; i++) {
            l.add(input[i]);
        }
        Collections.sort(l);
        for (int i = 0; i < l.size(); i++) {
            input[i] = l.get(i);
        }
    }

    public static void main(String[] args) {
        new AoColor();

    }

}
