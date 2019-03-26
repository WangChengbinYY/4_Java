package com.rs;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.math.Vector2D;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.simple.SimpleFeatureSource;
import org.opengis.feature.simple.SimpleFeature;

import java.io.*;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import static java.lang.Math.*;
import static java.lang.StrictMath.exp;

/**
 * Created by cwj
 * Date: 2018/8/10 0010
 * Modified by cwj
 * Date: 2018/8/10 0010
 * Abstract:
 */

public class rsConverter {
    public static final double M_PI = 3.14159265358979323846;

    /**
     * @param lat
     * @param lon
     * @return
     * 利用米勒投影计算，经纬度坐标转平面坐标
     */
    public static double[] millierConvert2XY(double lat, double lon){
        double L = 6381372 * M_PI * 2;//地球周长
        double W = L;// 平面展开后，x轴等于周长
        double H = L / 2;// y轴约等于周长一半
        double mill = 2.3;// 米勒投影中的一个常数，范围大约在正负2.3之间
        double x = lon * M_PI / 180;// 将经度从度数转换为弧度
        double y = lat * M_PI / 180;// 将纬度从度数转换为弧度
        y = 1.25 * log(tan(0.25 * M_PI + 0.4 * y));// 米勒投影的转换
        // 弧度转为实际距离
        x = (W / 2) + (W / (2 * M_PI)) * x;
        y = (H / 2) - (H / (2 * mill)) * y;
        double[] result = new double[2];
        result[0] = x;
        result[1] = y;
        return result;

    }

    //替代上面的    millierConvert2XY
    /**
     * @param lon 经度
     * @param lat 纬度
     * @param h   高程
     * @return  高斯坐标 X Y Z
     * 利用米勒投影计算，经纬度坐标转平面坐标
     */
    public static double[] LBHtoXYZ( double lon,double lat, double h){
        double e,f,a,b,N,mlat,mlon,X,Y,Z;
        a = 6378137.0;
        b = 6356752.3142451795;
        f = (a-b)/a;
        e = Math.sqrt(2*f-f*f);
        mlat = lat * M_PI / 180;
        mlon = lon * M_PI / 180;
        N = a/Math.sqrt(1-e*e*Math.sin(mlat)*Math.sin(mlat));

        X = (N+h)*Math.cos(mlat)*Math.cos(mlon);
        Y = (N+h)*Math.cos(mlat)*Math.sin(mlon);
        Z = (N*(1-e*e)+h)*Math.sin(mlat);
        double[] result = new double[3];
        result[0] = X;
        result[1] = Y;
        result[2] = Z;
        return result;
    }


    /**
     * @param x
     * @param y
     * @return
     * 利用米勒投影计算，平面坐标转经纬度坐标
     */
    public static double[] millierConvert2LonLat(double x, double y){
        double L = 6381372 * M_PI * 2;//地球周长
        double W = L;// 平面展开后，x轴等于周长
        double H = L / 2;// y轴约等于周长一半
        double mill = 2.3;// 米勒投影中的一个常数，范围大约在正负2.3之间
        double lat;
        lat = ((H / 2 - y) * 2 * mill) / (1.25 * H);
        lat = ((atan(exp(lat)) - 0.25 * M_PI) * 180) / (0.4 * M_PI);
        double lon;
        lon = (x - W / 2) * 360 / W;
        double[] result = new double[2];
        result[0] = lon;
        result[1] = lat;
        return result;
    }

    //替代上面的 millierConvert2LonLat
    /**
     * @param x
     * @param y
     * @param z
     * @return lon lat h
     * 利用米勒投影计算，平面坐标转经纬度坐标
     */
    public static double[] XYZtoLBH(double x, double y,double z){
        double e2,exing2,f,a,b,N,p,Cita,mlat,mlon,h;
        a = 6378137.0;
        b = 6356752.3142451795;
        f = (a-b)/a;
        e2 = 2*f-f*f;
        p = Math.sqrt(x*x+y*y);
        Cita = Math.atan((z*a)/(p*b));
        exing2 = (a*a - b*b)/(b*b);

        mlat = Math.atan((z+exing2*b*Math.sin(Cita)*Math.sin(Cita)*Math.sin(Cita))/(p-e2*a*Math.cos(Cita)*Math.cos(Cita)*Math.cos(Cita)));
        mlon = Math.atan2(y,x);
        N = a/Math.sqrt(1-e2*Math.sin(mlat)*Math.sin(mlat));
        h = p/Math.cos(mlat) - N;
        mlat = mlat * 180 / M_PI;
        mlon = mlon * 180 / M_PI;

        double[] result = new double[3];
        result[0] = mlon;
        result[1] = mlat;
        result[2] = h;
        return result;
    }


    /**
     *
     * @param long1 经度1
     * @param lat1 维度1
     * @param long2 经度2
     * @param lat2 纬度2
     * @return
     */
    public static double getDistance(double long1, double lat1,double h1, double long2, double lat2,double h2) {
//        double a, b, R;
//        R = 6378137; // 地球半径
//        lat1 = lat1 * Math.PI / 180.0;
//        lat2 = lat2 * Math.PI / 180.0;
//        a = lat1 - lat2;
//        b = (long1 - long2) * Math.PI / 180.0;
//        double d;
//        double sa2, sb2;
//        sa2 = Math.sin(a / 2.0);
//        sb2 = Math.sin(b / 2.0);
//        d = 2 * R * Math.asin(Math.sqrt(sa2 * sa2 + Math.cos(lat1) * Math.cos(lat2) * sb2 * sb2));
//        return d;
        double XYZ1[] = LBHtoXYZ(long1,lat1,h1);
        double XYZ2[] = LBHtoXYZ(long2,lat2,h2);
        double distance;
        distance = Math.sqrt((XYZ1[0]-XYZ2[0])*(XYZ1[0]-XYZ2[0]) +  (XYZ1[1]-XYZ2[1])*(XYZ1[1]-XYZ2[1]) + (XYZ1[2]-XYZ2[2])*(XYZ1[2]-XYZ2[2]));
        return distance;
    }

    /**
     *
     * @param lat_a 纬度1
     * @param lng_a 经度1
     * @param lat_b 纬度2
     * @param lng_b 经度2
     * @return
     */
    public static double getAngle(double lat_a, double lng_a, double lat_b, double lng_b) {
/*
        double lat_aH = lat_a*Math.PI/180.0;
        double lng_aH = lng_a*Math.PI/180.0;
        double lat_bH = lat_b*Math.PI/180.0;
        double lng_bH = lng_b*Math.PI/180.0;

        double x = Math.sin(lng_bH - lng_aH)*Math.cos(lat_bH);
        double y = cos(lat_aH)*sin(lat_bH) - sin(lat_aH)*cos(lat_bH)*cos(lng_bH-lng_aH);
        double A = Math.atan2(x,y);
        double Bearing = A - 360.0*((int)(A/360.0));
        return Math.toDegrees(Bearing);

*/


        double delta_x = (lng_b - lng_a)*Math.PI/180.0;
        double delta_y = (lat_b - lat_a)*Math.PI/180.0;
        //增加 经度差 *  cos(纬度)
        double mMeanLat = (lat_b + lat_a)*Math.PI/180.0/2.0;
        delta_x = delta_x * Math.cos(mMeanLat);

        double cita = 0.0;
        if(Math.abs(delta_x) < 0.00000000001) {
            if(delta_y > 0.0) {
                return 0.0;
            }
            if(delta_y < 0.0) {
                return 180.0;
            }
        }
        if(Math.abs(delta_y) < 0.0000000001) {
            if(delta_x > 0.0) {
                return 90.0;
            }
            if(delta_x < 0.0) {
                return 270.0;
            }
        }
        //此处 delta_y 肯定不等于 0.0
        double Cita = 0.0;
        Cita = Math.atan(Math.abs(delta_x/delta_y));
        Cita = Math.toDegrees(Cita);
        if((delta_x > 0.0)  &&  (delta_y > 0.0))
            return Cita;
        if((delta_x > 0.0)  &&  (delta_y < 0.0))
            return 180.0 - Cita;
        if((delta_x < 0.0)  &&  (delta_y < 0.0))
            return 180.0 + Cita;
        if((delta_x < 0.0)  &&  (delta_y > 0.0))
            return 360.0 - Cita;

        return 0.0;

    }

    /**
     * @param obj
     * @param lstObj
     * @return
     * 判断集合里面是否有相同点
     */
    public static boolean isSamePoint(HashMap<String, Object> obj, List<HashMap<String, Object>> lstObj){
        boolean bSame = false;
        Coordinate cd = (Coordinate)obj.get("coordinate");
        int cdCode = (int)obj.get("code");
        for (int i=0; i<lstObj.size(); i++){
            HashMap<String, Object> item = lstObj.get(i);
            Coordinate cdItem = (Coordinate)item.get("coordinate");
            int code = (int)item.get("code");
            if (cdItem.equals3D(cd) && (cdCode != code)){
                bSame = true;
                break;
            }
        }
        return  bSame;
    }

    /**
     * @param srcCoords
     * @param deltaDis
     * @return
     * 根据指定点间距，计算并插值点
     */
    public static Coordinate[] calcPosInterpo(Coordinate[] srcCoords, double deltaDis){
        int nCount = srcCoords.length;
        int i = 0;
        double lastDis = 0;
        ArrayList<Coordinate> allPos = new ArrayList<>();
//        while(i <= nCount-1){
//            Coordinate curPoint = srcCoords[i];
//            allPos.add(curPoint);
//            i++;
//        }



      //如果输入的点串只有一个点，直接返回!
        if(nCount == 1)
        {
            return  srcCoords;
        }

        //加入第一个点
        allPos.add(srcCoords[0]);
        //计算中间的点
        i = 1;
        while ( i <= nCount -1)
        {
            Coordinate curPoint = srcCoords[i-1];
            Coordinate fromPoint = srcCoords[i];
            double curDis = getDistance(curPoint.x, curPoint.y,curPoint.z, fromPoint.x, fromPoint.y,fromPoint.z);
            int mNumber = (int)(curDis/deltaDis);
            if(mNumber >= 1)
            {
                double XYZ_Start[] = LBHtoXYZ(curPoint.x, curPoint.y,curPoint.z);
                double XYZ_End[] = LBHtoXYZ(fromPoint.x, fromPoint.y,fromPoint.z);
                double[] XYZ_Dis = new double[3];
                XYZ_Dis[0] = (XYZ_End[0] - XYZ_Start[0]);
                XYZ_Dis[1] = (XYZ_End[1] - XYZ_Start[1]);
                XYZ_Dis[2] = (XYZ_End[2] - XYZ_Start[2]);
                int j=0;
                for(j=0;j<mNumber;j++) {
                    Coordinate middlePoint = new Coordinate(0.0,0.0,0.0);
                    double[] XYZ_Middle = new double[3];
                    XYZ_Middle[0] = XYZ_Start[0] + ((j+1) * deltaDis / curDis) * XYZ_Dis[0];
                    XYZ_Middle[1] = XYZ_Start[1] + ((j+1) * deltaDis / curDis) * XYZ_Dis[1];
                    XYZ_Middle[2] = XYZ_Start[2] + ((j+1) * deltaDis / curDis) * XYZ_Dis[2];
                    double[] LBH_Middle = XYZtoLBH(XYZ_Middle[0], XYZ_Middle[1], XYZ_Middle[2]);
                    middlePoint.x = LBH_Middle[0];
                    middlePoint.y = LBH_Middle[1];
                    middlePoint.z = LBH_Middle[2];
                    //加入中间插值点
                    allPos.add(middlePoint);

                    Coordinate tempFore = allPos.get(allPos.size()-2);
                    double tempDis = getDistance(tempFore.x,tempFore.y,tempFore.z,middlePoint.x,middlePoint.y,middlePoint.z);
                    if (tempDis > 2.0)
                    {
                        double fuck;
                        fuck = 8.9;
                    }

                }
            }
            allPos.add(fromPoint);
            i++;
        }

        Coordinate[] posArrary = new Coordinate[allPos.size()];
        posArrary = allPos.toArray(posArrary);
        return  posArrary;
    }

    public static void convertLink2BinaryPoints(String flPath, String flTargetPath, String Encoding, int regionNum,double mDistance){
        long startTime= System.currentTimeMillis();
        int num= 0;

        ShapefileDataStore dsShp = null;
        SimpleFeatureIterator ftIter= null;
        try {
            dsShp = (ShapefileDataStore) FileDataStoreFinder.getDataStore(new File(flPath));
            dsShp.setCharset(Charset.forName(Encoding));
            SimpleFeatureSource shp = dsShp.getFeatureSource();
            ftIter = shp.getFeatures().features();

            if (ftIter == null){
                System.out.println(String.format("打开[%s]文件出错", flPath));
            }
            System.out.println(String.format("打开[%s]文件成功", flPath));

            // 获取所有路网线段的首尾点，用于判断提取点是否是交叉点
            SimpleFeatureIterator ite = shp.getFeatures().features();
            List<HashMap<String, Object>> lstObj = getBeginEndPoints(ite);
            ite.close();

            File file = new File(flTargetPath);
            OutputStream out = new FileOutputStream(file);

            int code = 1;
            while (ftIter.hasNext()){
                SimpleFeature ft = ftIter.next();
                String id = ft.getID();
                Geometry shape = (Geometry) ft.getDefaultGeometry();
                if (shape.isEmpty()|| !(shape.isValid())){
                    System.out.println(String.format("[%s]记录图形为空有图形无效", id));
                    continue;
                }

                System.out.println(String.format("正在转换[%s]的记录......", id));
                /*int code = Integer.parseInt(ft.getAttribute("CODE").toString());*/

                // 对道路中心2米
               Coordinate[] interpoCoors = calcPosInterpo(shape.getCoordinates(), mDistance);
//                Coordinate[] interpoCoors = shape.getCoordinates();

                convert2NavPoint(out, interpoCoors, regionNum, code++, lstObj);
                num++;
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (ftIter != null){
                ftIter.close();
            }
            if (dsShp != null){
                dsShp.dispose();
            }
        }

        long endTime=System.currentTimeMillis();
        System.out.println(String.format("转换[%d]记录------耗时: [%d]s", num, (endTime-startTime)/1000));
    }

    public static List<HashMap<String, Object>> getBeginEndPoints(SimpleFeatureIterator ite){
        List<HashMap<String, Object>> lstPtObj = new ArrayList<>();
        try {
            while (ite.hasNext()){
                SimpleFeature ft = ite.next();
                Geometry shape = (Geometry) ft.getDefaultGeometry();
                if (shape.isEmpty()|| !(shape.isValid())){
                    continue;
                }

                int code = Integer.parseInt(ft.getAttribute("CODE").toString());
                Coordinate[] cds = shape.getCoordinates();
                HashMap<String, Object> itemB = new HashMap<>();
                itemB.put("coordinate", cds[0]);
                itemB.put("code", code);
                lstPtObj.add(itemB);

                HashMap<String, Object> itemE = new HashMap<>();
                itemE.put("coordinate", cds[cds.length -1]);
                itemE.put("code", code);
                lstPtObj.add(itemE);
            }
        }catch (Exception e){
            e.printStackTrace();
        }
        return lstPtObj;
    }

    /**
     * 解析Geometry中的坐标点,转成输出格式的路网点
     */
    public static void convert2NavPoint(OutputStream out, Coordinate[] cds, int regNum, int code, List<HashMap<String, Object>> lstObj){
        try {
            for (int i=0; i<cds.length; i++){
                Coordinate cd = cds[i];
                navPointRec rec = new navPointRec();
                rec.regionNum = regNum;
                rec.roadNum = code;
                rec.roadSideNum = 0x00;
                rec.roadPointNum = i+1;
                rec.crossFlag = 0x00;

                HashMap<String,Object> obj = new HashMap<>();
                obj.put("coordinate", cd);
                obj.put("code", code);
                if ((i == 0 || i== cds.length -1) && isSamePoint(obj, lstObj)){
                    rec.crossFlag = 0x01;
                }
                rec.borderFlag = 0x00;
                rec.borderRelNum = 0;
                rec.longitude = cd.x;
                rec.latitude = cd.y;
                rec.altitude = (float) cd.z;

                if (i == cds.length -1)
                {
                    // 最后一个点取倒数第二点计算的方向
                    rec.direction = (float) getAngle(cds[i-1].y, cds[i-1].x, cd.y, cd.x);
                }
                else
                {
                    if(i == 0)
                    {
                        if((Math.abs(cd.y - cds[i+1].y) < 0.0000000001) && (Math.abs(cd.x - cds[i+1].x) < 0.0000000001))
                        {
                            continue;
                        }
                        else
                        {
                            rec.direction = (float) getAngle(cd.y, cd.x, cds[i+1].y, cds[i+1].x);
                        }
                    }
                    else
                    {
                        if((Math.abs(cd.y - cds[i+1].y) < 0.0000000001) && (Math.abs(cd.x - cds[i+1].x) < 0.0000000001))
                        {
                            rec.direction =  (float) getAngle(cds[i-1].y, cds[i-1].x, cd.y, cd.x);
                        }
                        else
                        {
                            rec.direction = (float) getAngle(cd.y, cd.x, cds[i+1].y, cds[i+1].x);
                        }
                     }
                }

                byte[] recBuf = byteMergerAll(
                        LittleByteUtil.getIntBytes(rec.regionNum),
                        LittleByteUtil.getIntBytes(rec.roadNum),
                        LittleByteUtil.getbyteBytes(rec.roadSideNum),
                        LittleByteUtil.getIntBytes(rec.roadPointNum),
                        LittleByteUtil.getbyteBytes(rec.crossFlag),
                        LittleByteUtil.getbyteBytes(rec.borderFlag),
                        LittleByteUtil.getIntBytes(rec.borderRelNum),
                        LittleByteUtil.getDoubleBytes(rec.longitude),
                        LittleByteUtil.getDoubleBytes(rec.latitude),
                        LittleByteUtil.getFloatBytes(rec.altitude),
                        LittleByteUtil.getFloatBytes(rec.direction));
                out.write(recBuf);

                double mDistance = 0.0;
                if(i < cds.length -1) {
                    mDistance = getDistance( cds[i].x,cds[i].y,cds[i].z,cds[i+1].x, cds[i+1].y,cds[i+1].z);
                }
                String tempString = String.valueOf(mDistance);
                System.out.println(String.format("写入point[%d]: %s,%s", i, rec.toString(),tempString));
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void readBinaryPoints(String flPath){
        try {
            int recByteSize = 43;
            byte[] buf = getBytes4File(flPath);
            int count = buf.length/recByteSize;
            for (int i=0; i<count; i++){
                byte[] rec = getSubBytes(buf, i*recByteSize, recByteSize);
                navPointRec navPt = new navPointRec();
                navPt.regionNum = LittleByteUtil.getInt(getSubBytes(rec, 0, 4));
                navPt.roadNum = LittleByteUtil.getInt(getSubBytes(rec, 4, 4));
                navPt.roadSideNum = LittleByteUtil.getbyte(getSubBytes(rec, 8, 1));
                navPt.roadPointNum = LittleByteUtil.getInt(getSubBytes(rec, 9, 4));
                navPt.crossFlag = LittleByteUtil.getbyte(getSubBytes(rec, 13, 1));
                navPt.borderFlag = LittleByteUtil.getbyte(getSubBytes(rec, 14, 1));
                navPt.borderRelNum = LittleByteUtil.getInt(getSubBytes(rec, 15, 4));
                navPt.longitude = LittleByteUtil.getDouble(getSubBytes(rec, 19, 8));
                navPt.latitude = LittleByteUtil.getDouble(getSubBytes(rec, 27, 8));
                navPt.altitude = LittleByteUtil.getFloat(getSubBytes(rec, 35, 4));
                navPt.direction = LittleByteUtil.getFloat(getSubBytes(rec, 39, 4));

                System.out.println(String.format("读取point[%d]: %s", i, navPt.toString()));
            }
        }catch (IOException e){
            e.printStackTrace();
        }
    }

    /**
     * 读取二进制文件并且写入数组里
     * @param filePath
     * @return
     * @throws IOException
     * @throws FileNotFoundException
     */
    public static byte[] getBytes4File(String filePath) throws IOException {

        InputStream in = null;
        BufferedInputStream buffer = null;
        DataInputStream dataIn = null;
        ByteArrayOutputStream bos = null;
        DataOutputStream dos = null;
        byte[] bArray = null;
        try {
            in = new FileInputStream(filePath);
            buffer = new BufferedInputStream(in);
            dataIn = new DataInputStream(buffer);
            bos = new ByteArrayOutputStream();
            dos = new DataOutputStream(bos);
            byte[] buf = new byte[1024];
            while (true) {
                int len = dataIn.read(buf);
                if (len < 0)
                    break;
                dos.write(buf, 0, len);
            }
            bArray = bos.toByteArray();

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        } finally {

            if (in != null)
                in.close();
            if (dataIn != null)
                dataIn.close();
            if (buffer != null)
                buffer.close();
            if (bos != null)
                bos.close();
            if (dos != null)
                dos.close();
        }

        return bArray;
    }

    /**
     * 从一个byte[]数组中截取一部分
     * @param src
     * @param begin
     * @param count
     * @return
     */
    public static byte[] getSubBytes(byte[] src, int begin, int count) {
        byte[] bs = new byte[count];
        for (int i=begin;i<begin+count; i++) bs[i-begin] = src[i];
        return bs;
    }

    /***
     * 合并字节数组
     *
     * @param a
     * @return
     */
    public static byte[] byteMergerAll(byte[]... a) {
        // 合并完之后数组的总长度
        int index = 0;
        int sum = 0;
        for (int i = 0; i < a.length; i++) {
            sum = sum + a[i].length;
        }
        byte[] result = new byte[sum];
        for (int i = 0; i < a.length; i++) {
            int lengthOne = a[i].length;
            if (lengthOne == 0) {
                continue;
            }
            // 拷贝数组
            System.arraycopy(a[i], 0, result, index, lengthOne);
            index = index + lengthOne;
        }
        return result;
    }
}
