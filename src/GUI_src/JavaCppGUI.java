import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;

public class JavaCppGUI{
    private JFrame frame;
    private JTextField inputField;
    private JTextArea outputArea;
    private JButton runButton;

    public JavaCppGUI(){
        frame = new JFrame("Java GUI for C++ Program");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(1600, 900);
        frame.setLayout(new BorderLayout());

        inputField = new JTextField();
        inputField.setEditable(true);
        // inputField.setFont();
        outputArea = new JTextArea();
        outputArea.setEditable(false);
        runButton = new JButton("Run C++ Program");

        runButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                runCppProgram();
            }
        });

        frame.add(inputField, BorderLayout.NORTH);
        frame.add(new JScrollPane(outputArea), BorderLayout.CENTER);
        frame.add(runButton, BorderLayout.SOUTH);
        frame.setVisible(true);
    }

    private void runCppProgram(){
        try {
            FileWriter writer = new FileWriter("./input.txt");
            writer.write(inputField.getText());
            writer.close();

            Process process = Runtime.getRuntime().exec("./cpp_program");

            process.waitFor();

            BufferedReader reader = new BufferedReader(new FileReader("./output.txt"));
            String line;
            StringBuilder output = new StringBuilder();
            while( (line = reader.readLine()) != null ){
                output.append(line).append("\n");
            }
            reader.close();

            outputArea.setText(output.toString());
        } catch(Exception ex){
            ex.printStackTrace();
            outputArea.setText("Error running C++ program.");
        }
    }

    public static void main(String[] args){
        SwingUtilities.invokeLater(JavaCppGUI::new);
    }
}
