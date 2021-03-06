package com.rs;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;
import javax.swing.filechooser.FileFilter;
import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.IOException;

class Visualize extends JFrame implements ActionListener {

    private Font titleFont= new Font("Microsoft YaHei", Font.BOLD, 14);
    private Font textFont= new Font("Microsoft YaHei", Font.BOLD, 12);
    private EmptyBorder border= new EmptyBorder(0,0,0,0);

    private JTextField inputFile;
    private JButton inputButton;
    private JTextField outputPath;
    private JButton outputButton;
    private JComboBox<String> encoding;
    private JTextField roadNetwork;
    private JButton execute;
    private ConsoleTextArea consoleText;
    private JTextField interpositionPoint;

    Visualize(){
        super();
        setFont(titleFont);
        setTitle("道路路网提取");
        double lx = Toolkit.getDefaultToolkit().getScreenSize().getWidth();
        double ly = Toolkit.getDefaultToolkit().getScreenSize().getHeight();
        setBounds((int)(lx /2- 160), (int)(ly /2-225), 430, 500);
        setLayout(null);

        add(filePathText());
        add(filePath());
        add(fileButton());
        add(outputPathText());
        add(outputPath());
        add(pathButton());
        add(roadNetworkText());
        add(interpostionPointText());
        add(encoding());
        add(executeButton());
        add(roadNetwork());
        add(interpositionPoint());
        add(encodingText());

        try {
            add(getConsole());
        } catch (IOException e) {
            e.printStackTrace();
        }

        setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        setVisible(true);
    }

    private JTextField filePathText(){
        JTextField inputFileText= new JTextField();
        inputFileText.setText("导入文件：");
        inputFileText.setLayout(null);
        inputFileText.setEditable(false);
        inputFileText.setBorder(border);
        inputFileText.setFont(titleFont);
        inputFileText.setBounds(15, 10, 90, 20);
        return inputFileText;
    }

    private JTextField filePath(){
        inputFile= new JTextField();
        inputFile.setLayout(null);
        inputFile.setEditable(false);
        inputFile.addActionListener(this);
        inputFile.setBackground(Color.WHITE);
        inputFile.setFont(textFont);
        inputFile.setBounds(15, 30, 265, 25);
        return inputFile;
    }

    private JButton fileButton(){
        inputButton= new JButton();
        inputButton.setText("打开文件");
        inputButton.setFont(textFont);
        inputButton.setLayout(null);
        inputButton.addActionListener(this);
        inputButton.setBounds(290, 30, 90, 25);
        return inputButton;
    }

    private JTextField outputPathText(){
        JTextField outputPathText= new JTextField();
        outputPathText.setText("输出路径：");
        outputPathText.setLayout(null);
        outputPathText.setEditable(false);
        outputPathText.setBorder(border);
        outputPathText.setFont(titleFont);
        outputPathText.setBounds(15, 65, 90, 20);
        return outputPathText;
    }

    private JTextField outputPath(){
        outputPath= new JTextField();
        outputPath.setLayout(null);
        outputPath.setEditable(false);
        outputPath.addActionListener(this);
        outputPath.setBackground(Color.WHITE);
        outputPath.setFont(textFont);
        outputPath.setBounds(15, 85, 265, 25);
        return outputPath;
    }

    private JButton pathButton(){
        outputButton= new JButton();
        outputButton.setText("打开目录");
        outputButton.setFont(textFont);
        outputButton.setLayout(null);
        outputButton.addActionListener(this);
        outputButton.setBounds(290, 85, 90, 25);
        return outputButton;
    }

    private JTextField encodingText(){
        JTextField encodingText= new JTextField();
        encodingText.setText("导入文件编码：");
        encodingText.setLayout(null);
        encodingText.setEditable(false);
        encodingText.setBorder(border);
        encodingText.setFont(titleFont);
        encodingText.setBounds(15, 130, 90, 20);
        return encodingText;
    }

    private JComboBox encoding(){
        encoding= new JComboBox<>();
        encoding.setFont(textFont);
        encoding.setBackground(Color.WHITE);
        encoding.setBounds(115, 130, 70, 25);
        encoding.addItem("UTF-8");
        encoding.addItem("GBK");
        encoding.setSelectedItem("UTF8");
        return encoding;
    }

    private JTextField roadNetworkText(){
        JTextField roadNetworkText= new JTextField();
        roadNetworkText.setText("道路路网编号：");
        roadNetworkText.setLayout(null);
        roadNetworkText.setEditable(false);
        roadNetworkText.setBorder(border);
        roadNetworkText.setFont(titleFont);
        roadNetworkText.setBounds(200, 130, 90, 25);
        return roadNetworkText;
    }

    private JTextField interpostionPointText(){
        JTextField interpostionPointText= new JTextField();
        interpostionPointText.setText("道路插点间隔(>1.0m)：");
        interpostionPointText.setLayout(null);
        interpostionPointText.setEditable(false);
        interpostionPointText.setBorder(border);
        interpostionPointText.setFont(titleFont);
        interpostionPointText.setBounds(15, 180, 150, 25);
        return interpostionPointText;
    }

    private JTextField interpositionPoint(){
        interpositionPoint= new JTextField();
        interpositionPoint.setLayout(null);
        interpositionPoint.setFont(textFont);
        interpositionPoint.setBounds(170, 180, 70, 25);
        interpositionPoint.addKeyListener(new KeyAdapter(){
            public void keyTyped(KeyEvent e){
                int keyChar = e.getKeyChar();
                if((!(keyChar >= KeyEvent.VK_0 && keyChar <= KeyEvent.VK_9)) && (keyChar !=KeyEvent.VK_PERIOD))
                    e.consume();    //关键，屏蔽掉非法输入
            }
        }); //限制只能输入数字 和 小数点
        return interpositionPoint;
    }




    private JTextField roadNetwork(){
        roadNetwork= new JTextField();
        roadNetwork.setLayout(null);
        roadNetwork.setFont(textFont);
        roadNetwork.setBounds(300, 130, 70, 25);
        roadNetwork.addKeyListener(new KeyAdapter(){
            public void keyTyped(KeyEvent e){
                int keyChar = e.getKeyChar();
                if(!(keyChar >= KeyEvent.VK_0 && keyChar <= KeyEvent.VK_9))
                    e.consume();    //关键，屏蔽掉非法输入
            }
        }); //限制只能输入数字
        return roadNetwork;
    }

    private JTextArea getConsole() throws IOException {
        PopupMenu pMenu=new PopupMenu();
        MenuItem mItemClear=new MenuItem("clear");
        MouseAdapter mouseAdapter=new MouseAdapter(){
            public void mouseClicked(MouseEvent event)
            {
                if(event.getButton()==MouseEvent.BUTTON3)
                {
                    pMenu.show(consoleText,event.getX(),event.getY());
                }
            }
        };
        ActionListener menuAction=new ActionListener() {
            public void actionPerformed(ActionEvent e)
            {
                MenuItem item=(MenuItem)e.getSource();
                if(item==mItemClear) {
                    consoleText.setText("");
                }
            }
        };

        consoleText= new ConsoleTextArea();
        consoleText.setEditable(false);
        consoleText.add(pMenu);
        consoleText.addMouseListener(mouseAdapter);
        pMenu.add(mItemClear);
        mItemClear.addActionListener(menuAction);
        JScrollPane jsp= new JScrollPane(consoleText);
        jsp.setBounds(0, 0, 380, 185);
        jsp.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        jsp.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        JTextArea console= new JTextArea();
        console.setEditable(false);
        console.setFont(textFont);
        console.setBounds(10, 220, 380, 185);
        console.setBorder(new LineBorder(Color.lightGray));
        console.add(jsp);
        return console;
    }

    private JButton executeButton(){
        execute= new JButton();
        execute.setText("确定");
        execute.setLayout(null);
        execute.setFont(titleFont);
        execute.setBounds(170, 420, 80, 25);
        execute.addActionListener(this);
        return execute;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        if (e.getSource().equals(inputButton)){
            JFileChooser fileSelect = new JFileChooser();
            fileSelect.setFileSelectionMode(JFileChooser.FILES_ONLY);
            fileSelect.setFileFilter(new FileFilter() {
                @Override
                public boolean accept(File f) {
                    return f.getName().endsWith(".shp") || f.isDirectory();
                }

                @Override
                public String getDescription() {
                    return "shape文件(*.shp)";
                }
            });
            int state = fileSelect.showOpenDialog(null);
            if (state!= 1) {
                File f = fileSelect.getSelectedFile();
                inputFile.setText(f.getAbsolutePath());
            }
        }

        if (e.getSource().equals(outputButton)){
            JFileChooser fileSelect = new JFileChooser();
            fileSelect.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            int state = fileSelect.showOpenDialog(null);
            if (state!= 1) {
                File f = fileSelect.getSelectedFile();
                outputPath.setText(f.getAbsolutePath());
            }
        }

        if (e.getSource().equals(execute)){
            String input= inputFile.getText();
            assert input != null: "导入文件不能为空";
            String inputName= new File(input).getName();
            String outputName= inputName.substring(0, inputName.lastIndexOf("."))+ ".dat";
            String outputDir= outputPath.getText();
            assert outputDir != null: "导出路径不能为空";
            String output= outputDir+ File.separator+ outputName;
            String encode= String.valueOf(encoding.getSelectedItem());
            String network= roadNetwork.getText();
            assert network != null: "道路路网编号不能为空";

            String interpostionDistance= interpositionPoint.getText();
            assert network != null: "插值距离不能为空(>1.0m)";

            Float tempDistanc = Float.parseFloat(interpostionDistance);
            if(tempDistanc < 1.0){
                rsConverter.convertLink2BinaryPoints(input, output, encode, Integer.parseInt(network),1.0);
            }else{
                rsConverter.convertLink2BinaryPoints(input, output, encode, Integer.parseInt(network),Float.parseFloat(interpostionDistance));
            }



        }

    }

}
