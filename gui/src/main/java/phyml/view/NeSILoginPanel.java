package phyml.view;

import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;
import grith.gsindl.SLCS;
import grith.sibboleth.ShibLoginPanel;

import javax.swing.*;

/**
 * Project: grisu
 * <p/>
 * Written by: Markus Binsteiner
 * Date: 19/08/13
 * Time: 2:57 PM
 */
public class NeSILoginPanel {
    private JPanel panel1;
    private ShibLoginPanel shibLoginPanel1;
    private JButton loginButton;


    public JPanel getPanel() {
        return panel1;
    }

    private void createUIComponents() {
        shibLoginPanel1 = new ShibLoginPanel(SLCS.DEFAULT_SLCS_URL);
    }

    {
// GUI initializer generated by IntelliJ IDEA GUI Designer
// >>> IMPORTANT!! <<<
// DO NOT EDIT OR ADD ANY CODE HERE!
        $$$setupUI$$$();
    }

    /**
     * Method generated by IntelliJ IDEA GUI Designer
     * >>> IMPORTANT!! <<<
     * DO NOT edit this method OR call it in your code!
     *
     * @noinspection ALL
     */
    private void $$$setupUI$$$() {
        createUIComponents();
        panel1 = new JPanel();
        panel1.setLayout(new FormLayout("fill:d:grow", "center:d:grow,top:4dlu:noGrow,center:d:noGrow,top:4dlu:noGrow,top:d:grow,top:4dlu:noGrow,center:max(d;4px):noGrow"));
        final JLabel label1 = new JLabel();
        label1.setText("You are not logged in.");
        CellConstraints cc = new CellConstraints();
        panel1.add(label1, cc.xy(1, 3));
        panel1.add(shibLoginPanel1, cc.xy(1, 5, CellConstraints.FILL, CellConstraints.FILL));
        loginButton = new JButton();
        loginButton.setHorizontalAlignment(4);
        loginButton.setHorizontalTextPosition(0);
        loginButton.setText("Login");
        panel1.add(loginButton, cc.xy(1, 7, CellConstraints.RIGHT, CellConstraints.DEFAULT));
        final JSeparator separator1 = new JSeparator();
        panel1.add(separator1, cc.xy(1, 1, CellConstraints.FILL, CellConstraints.FILL));
    }

    /**
     * @noinspection ALL
     */
    public JComponent $$$getRootComponent$$$() {
        return panel1;
    }
}