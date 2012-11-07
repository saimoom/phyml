
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Rectangle;
import javax.swing.JPanel;

/**
 * Extends a JPanel to provide a possibility of a Separator with a different
 * colour or a label.
 * 
 * @author Christoph Knapp
 * 
 */
public class Separator extends JPanel {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private String text;
	private Color forground;
	boolean isEmpty;

	/**
	 * Constructor method to provide a possibility of a Separator with a
	 * different colour or a label.
	 * 
	 * @param txt
	 *            String : sets the text of the label
	 * @param background
	 *            Color : sets the background Color
	 * @param forground
	 *            Color : sets the forground Color
	 */
	public Separator(String txt, Color background, Color forground) {
		isEmpty = false;
		text = txt;
		if (background != null) {
			setBackground(background);
		}
		this.forground = forground;
	}

	/**
	 * Constructor method to provide a possibility of a Separator with a
	 * different colour or a label.
	 */
	public Separator() {
		isEmpty = false;
		text = "";
		this.forground = null;
	}

	/**
	 * Constructor method to provide a possibility of a Separator with a
	 * different colour or a label.
	 * 
	 * @param background
	 *            Color : sets the background Color
	 * @param forground
	 *            Color : sets the forground Color
	 */
	public Separator(Color background, Color forground) {
		isEmpty = false;
		text = "";
		if (background != null) {
			setBackground(background);
		}
		this.forground = forground;
	}

	/**
	 * Constructor method to provide a possibility of a Separator with a
	 * different colour or a label.
	 * 
	 * @param isEmpty
	 *            boolean : true if there is not a horizontal line, false if
	 *            there is.
	 */
	public Separator(boolean isEmpty) {
		this.isEmpty = isEmpty;
		text = "";
		this.forground = null;
	}

	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		if (!isEmpty) {
			if (forground != null) {
				g.setColor(forground);
			}
			if (!text.equals("")) {
				Rectangle r = g.getFontMetrics().getStringBounds(text, g)
						.getBounds();
				g.drawLine(5, this.getHeight() / 2,
						(int) (this.getWidth() / 2 - r.getWidth() / 2) - 5,
						this.getHeight() / 2);
				g.drawLine((int) (this.getWidth() / 2 + r.getWidth() / 2) + 5,
						this.getHeight() / 2, (int) (this.getWidth()) - 5,
						this.getHeight() / 2);
				g.drawString(text,
						(int) (this.getWidth() / 2 - r.getWidth() / 2),
						(int) (this.getHeight() / 2 + r.getHeight() / 2));
			} else {
				g.drawLine(5, this.getHeight() / 2, this.getWidth() - 5,
						this.getHeight() / 2);
			}
		}
	}
}
