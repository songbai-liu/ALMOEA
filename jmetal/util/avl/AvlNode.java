package jmetal.util.avl;

/**
 * Created with IntelliJ IDEA. User: Antonio J. Nebro Date: 08/07/13 Time: 15:46
 * To change this template use File | Settings | File Templates.
 */
public class AvlNode<T> {
	private AvlNode<T> left_;
	private AvlNode<T> right_;
	private AvlNode<T> parent_;

	private int height_;

	private AvlNode<T> closestNode_;

	private T item_;

	/**
	 * Constructor
	 * 
	 * @param item_
	 */
	public AvlNode(T item_) {
		this.left_ = null;
		this.right_ = null;
		this.parent_ = null;
		height_ = 0;
		closestNode_ = null;

		this.item_ = item_;
	}

	public AvlNode getLeft() {
		return left_;
	}

	public void setLeft(AvlNode left) {
		this.left_ = left;
	}

	public AvlNode getParent() {
		return parent_;
	}

	public void setParent(AvlNode parent) {
		this.parent_ = parent;
	}

	public AvlNode getRight() {
		return right_;
	}

	public void setRight(AvlNode right) {
		this.right_ = right;
	}

	public T getItem() {
		return item_;
	}

	public void setItem(T item) {
		this.item_ = item;
	}

	public int getHeight() {
		return height_;
	}

	public void setHeight(int height) {
		this.height_ = height;
	}

	public void updateHeight() {
		if (!hasLeft() && !hasRight()) {
			height_ = 0;
		} else if (!hasRight()) {
			height_ = 1 + getLeft().getHeight();
		} else if (!hasLeft()) {
			height_ = 1 + getRight().getHeight();
		} else {
			height_ = 1 + Math.max(getLeft().getHeight(), getRight()
					.getHeight());
		}
	}

	public AvlNode<T> getClosestNode() {
		return closestNode_;
	}

	public void setClosestNode_(AvlNode<T> closestNode) {
		this.closestNode_ = closestNode;
	}

	public boolean hasParent() {
		return parent_ != null;
	}

	public boolean hasLeft() {
		return left_ != null;
	}

	public boolean hasRight() {
		return right_ != null;
	}

	public boolean isLeaf() {
		boolean result;
		if (!hasLeft() && !hasRight()) {
			result = true;
		} else
			result = false;

		return result;
	}

	public boolean hasOnlyALeftChild() {
		boolean result;
		if (hasLeft() && !hasRight()) {
			result = true;
		} else
			result = false;

		return result;
	}

	public boolean hasOnlyARightChild() {
		boolean result;
		if (hasRight() && !hasLeft()) {
			result = true;
		} else
			result = false;

		return result;
	}
}
