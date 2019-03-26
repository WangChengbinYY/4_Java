package com.rs;

/**
 * Created by cwj
 * Date: 2018/8/10 0010
 * Modified by cwj
 * Date: 2018/8/10 0010
 * Abstract:
 */

public class navPointRec {
    public int regionNum = 1;               // ·��������     4byte
    public int roadNum = 1;                 // ��·���        4byte
    public byte roadSideNum = 0x00;         // ��·����      1byte
    public int  roadPointNum = 1;           // ��·����      4byte
    public byte crossFlag = 0x00;           // ���������      1byte
    public byte borderFlag = 0x00;          // �߽������      1byte
    public int borderRelNum = 0;            // �߽��������    4byte
    public double longitude = 0.0;          // ����           8byte
    public double latitude = 0.0;           // γ��           8byte
    public float altitude = 0.0f;           // �߳�           4byte
    public float direction = 0.0f;          // ������Ϣ        4byte


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
