package com.rs;

import javax.swing.*;
import javax.swing.text.Document;
import java.io.*;

public class ConsoleTextArea extends JTextArea{

    public ConsoleTextArea(InputStream[] inStreams){
        for(int i = 0; i < inStreams.length; ++i)
            startConsoleReaderThread(inStreams[i]);
    }

    public ConsoleTextArea() throws IOException {
        final LoopedStreams ls = new LoopedStreams();
        // 重定向System.out和System.err
        PrintStream ps = new PrintStream(ls.getOutputStream());
        System.setOut(ps);
        System.setErr(ps);
        startConsoleReaderThread(ls.getInputStream());
    }

    private void startConsoleReaderThread(InputStream inStream) {
        final BufferedReader bufferedReader=
                new BufferedReader(new InputStreamReader(inStream));

        new Thread(new Runnable(){
            public void run() {
                StringBuffer stringBuffer=new StringBuffer();
                String s;
                Document document=getDocument();
                try {
                    while ((s=bufferedReader.readLine())!=null) {
                        boolean caretAtEnd=false;
                        caretAtEnd=getCaretPosition()==document.getLength()?true:false;
                        stringBuffer.setLength(0);
                        append(stringBuffer.append(s).append("\n").toString());
                        if (caretAtEnd) {
                            setCaretPosition(document.getLength());
                        }
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }).start();
    }
}
