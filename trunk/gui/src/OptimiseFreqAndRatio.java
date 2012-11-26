
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;

/**
 * This JPanel implements the components to specify whether the "Equilibrium
 * Frequency" is optimised or not and whether the "Transition/Transversion Ratio"
 * is estimated or user defined.
 * 
 * @author Christoph Knapp
 * @date 02/07/12
 */
public class OptimiseFreqAndRatio extends JPanel implements ActionListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
//	private JLabel optimiseLab;
//	private JComboBox optimiseBox;
	private JLabel ratioLab;
	private JComboBox<String> ratioBox;
	private CustomTextField ratioField;
	private final static String DNA = "DNA";
	private final static String AA = "AA";
	private boolean isDNA;
//	private String[] optimiseAray1;
	private String moleculeType;

	/**
	 * This JPanel implements the components to specify whether the "Equilibrium
	 * Frequency" is optimised or not and whether the "Transition/Transversion"
	 * ratio is estimated or user defined.
	 * 
	 * @param moleculeType
	 * String : Which molecule type the data is (DNA or AA).
	 */
	public OptimiseFreqAndRatio(String moleculeType) {
		this.moleculeType = moleculeType;
//		optimiseAray1 = new String[] { "yes", "no" };
		ratioLab = new JLabel("Transition/transversion ratio");
		ratioBox = new JComboBox<String>(new String[] { "estimated", "fixed" });
		ratioField = new CustomTextField("4.0");
		ratioField.setEnabled(false);
		if (moleculeType.equals(DNA)) {
			isDNA = true;
//			optimiseLab = new JLabel("Optimise Equ. Freq.");
//			optimiseLab.setToolTipText("Optimise equilibrium frequencies.");
//			optimiseBox = new JComboBox(optimiseAray1);
			isDNA = true;
		} else if (moleculeType.equals(AA)) {
			isDNA = false;
//			optimiseLab.setVisible(false);
			ratioLab.setVisible(false);
			ratioBox.setVisible(false);
			ratioField.setVisible(false);
		}
		ratioBox.addActionListener(this);
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.04);
		add(new JPanel());
		layout.setDimensions(0.01, 0.96);
		add(new JPanel());
		layout.setDimensions(0.98, 0.96);
		JPanel p2 = new JPanel();
		add(p2);
		CustomGridLayout l2 = new CustomGridLayout();
		p2.setLayout(l2);
		l2.setDimensions(0.691, 1);
		p2.add(ratioLab);
		l2.setDimensions(0.13, 1);
		p2.add(ratioBox);
		l2.setDimensions(0.05, 1);
		p2.add(new JPanel());
		l2.setDimensions(0.1, 1);
		p2.add(ratioField);
		layout.setDimensions(0.01, 0.96);
		add(new JPanel());
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == ratioBox) {
			if (ratioBox.getSelectedItem() == "fixed") {
				ratioField.setEnabled(true);
			} else {
				ratioField.setEnabled(false);
			}
		}
	}

	/**
	 * When the Moleculetype is changed in a DataType object the moleculetype
	 * needs to change here as well.
	 * 
	 * @param moleType
	 * String : Whether "DNA" or "AA" is selected in the DataType JPanel
	 */
	public void setMoleculeType(String moleType) {
		moleculeType = moleType;
		if (moleType.equals(DNA)) {
			isDNA = true;
//			optimiseLab.setVisible(true);
//			optimiseBox.setVisible(true);
			ratioLab.setVisible(true);
			ratioBox.setVisible(true);
			ratioField.setVisible(true);
		} else if (moleType.equals(AA)) {
			isDNA = false;
//			optimiseLab.setVisible(false);
//			optimiseBox.setVisible(false);
			ratioLab.setVisible(false);
			ratioBox.setVisible(false);
			ratioField.setVisible(false);
		}
	}

	/**
	 * Retrieves whether "DNA" or "AA" is set as molecule type.
	 * 
	 * @return
	 * boolean : returns true if molecule type is "DNA" and false when "AA".
	 */
	public boolean isDNA() {
		return isDNA;
	}

	/**
	 * Private class implementing a CustomTextField (extends JTextField) so that
	 * only double values can be entered.
	 */
	private class CustomTextField extends JTextField {
		/**
		 * default id
		 */
		private static final long serialVersionUID = 1L;

		/**
		 * Constructor only accepting Strings which can be converted into
		 * Double.
		 * 
		 * @param string
		 * String : Needs to be convertable into a Double value or it will
		 * not be accepted. Otherwise it sets the displayed default
		 * value of the TextField.
		 */
		public CustomTextField(String string) {
			super(string);
			if (!isNumber(string)) {
				this.setText("");
			}
		}

		protected Document createDefaultModel() {
			return new doubleOnlyDocument();
		}

		/**
		 * tests whether the input value can be converted into a Double value.
		 * 
		 * @param str
		 * String : TextField input String
		 * @return 
		 * boolean : true if str can be converted, false otherwise
		 */
		private boolean isNumber(String str) {
			try {
				Double.parseDouble(str);
				return true;
			} catch (NumberFormatException e) {
				if (str.equals(".") && !getTsTvRatio().contains(".")) {
					return true;
				}
				return false;
			}
		}

		/**
		 * Private class extending the PlainDocument document making sure only
		 * Double values can be typed into the JTextField
		 */
		private class doubleOnlyDocument extends PlainDocument {
			/**
			 * default id
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public void insertString(int offs, String str, AttributeSet a)
					throws BadLocationException {
				if (str == null || !isNumber(str)) {
					return;
				}
				super.insertString(offs, str, a);
			}
		}
	}

	/**
	 * Returns whether the TS/TV ratio is estimated or fixed.
	 * 
	 * @return 
	 * String : "e" for estimated or Double value from CustomTextField for fixed
	 * or "" otherwise.
	 */
	public String getTsTvRatio() {
		if (moleculeType.equals(DNA)) {
			if (ratioBox.getSelectedItem().toString().equals("estimated")) {
				return "e";
			}
			return ratioField.getText();
		}
		return "";
	}
	/**
	 * Enables or disables the Transition/Transversion Ratio dropdownmenu.
	 * 
	 * @param off
	 * boolean : if true disabled otherwise enabled;
	 */
	public void setRatioBoxOff(boolean off){
		ratioBox.setEnabled(!off);
	}
	/**
	 * Sets the index of ratioBox to 0.
	 */
	public void setRatioBoxDefault() {
		ratioBox.setSelectedIndex(0);
	}

}
