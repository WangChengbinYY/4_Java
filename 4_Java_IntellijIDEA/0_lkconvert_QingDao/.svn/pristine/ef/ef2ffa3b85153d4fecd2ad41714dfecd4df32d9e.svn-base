package com.rs;

/**
 * Created by cwj
 * Date: 2018/8/10 0010
 * Modified by cwj
 * Date: 2018/8/10 0010
 * Abstract:
 */

public class navPointRec {
    public int regionNum = 1;               // 路网区域编号     4byte
    public int roadNum = 1;                 // 道路编号        4byte
    public byte roadSideNum = 0x00;         // 道路侧编号      1byte
    public int  roadPointNum = 1;           // 道路点编号      4byte
    public byte crossFlag = 0x00;           // 交叉点属性      1byte
    public byte borderFlag = 0x00;          // 边界点属性      1byte
    public int borderRelNum = 0;            // 边界关联属性    4byte
    public double longitude = 0.0;          // 经度           8byte
    public double latitude = 0.0;           // 纬度           8byte
    public float altitude = 0.0f;           // 高程           4byte
    public float direction = 0.0f;          // 走向信息        4byte


    public String toString(){
        return String.format("%d, %d, %x, %d, %x, %x, %d, %f, %f, %f, %f",
                regionNum,
                roadNum,
                (int)roadSideNum,
                roadPointNum,
                (int)crossFlag,
                (int)borderFlag,
                borderRelNum,
                longitude,
                latitude,
                altitude,
                direction);
    }
}
