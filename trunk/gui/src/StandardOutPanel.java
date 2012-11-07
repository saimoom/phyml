
import javax.swing.DropMode;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;

/**
 * Extends JPanel and implements the components to display the standard output
 * of phyml.
 * 
 * @author Christoph Knapp
 * 
 */

public class StandardOutPanel extends JPanel {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private static JTextPane editorPane;

	/**
	 * Constructor method implements the components to display the standard
	 * output of phyml.
	 */
	public StandardOutPanel() {
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 1);
		editorPane = new JTextPane();
		editorPane.setEditable(false);
		editorPane.setDropMode(DropMode.INSERT);
		JScrollPane editorScrollPane = new JScrollPane(editorPane);
		editorScrollPane
				.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
		editorScrollPane
				.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		add(editorScrollPane);
	}

	/**
	 * Adds a new line to the standard output panel.
	 * 
	 * @param line
	 *            String : A line to add to the output panel
	 */
	public static void setInput(String line) {
		if(line.contains("\t")){
			System.out.println(line);
		}
		editorPane.setText(editorPane.getText() + line + "\n");
	}
}
